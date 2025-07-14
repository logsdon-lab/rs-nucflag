use coitrees::{GenericInterval, Interval};
use eyre::Context;
use itertools::Itertools;
use noodles::{
    bam, bgzf,
    core::{Position, Region},
    cram,
    sam::{
        alignment::record::{cigar::op::Kind, Cigar},
        Header,
    },
};
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::{ffi::OsStr, fmt::Debug, fs::File, path::Path};

use crate::config::Config;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum PileupMAPQFn {
    Median,
    #[default]
    Mean,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PileupInfo {
    pub n_cov: u32,
    pub n_mismatch: u32,
    pub n_indel: u32,
    pub n_softclip: u32,
    pub mapq: Vec<u8>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct PileupSummary {
    pub region: Region,
    pub pileups: Vec<PileupInfo>,
}

pub enum AlignmentFile {
    Cram(cram::io::IndexedReader<File>),
    Bam(bam::io::IndexedReader<bgzf::Reader<File>>),
}

impl PileupInfo {
    pub fn median_mapq(&self) -> Option<u8> {
        let length = self.mapq.len();
        let midpt = length / 2;
        if length % 2 == 0 {
            let midpt = midpt.checked_sub(1).map(|midpt| midpt..=midpt)?;
            Some(self.mapq.iter().sorted().get(midpt).sum::<u8>().div_ceil(2))
        } else {
            self.mapq.iter().sorted().nth(self.mapq.len() / 2).cloned()
        }
    }
    pub fn mean_mapq(&self) -> eyre::Result<u8> {
        let Some(length) = TryInto::<u32>::try_into(self.mapq.len())
            .ok()
            .filter(|l| *l > 0)
        else {
            return Ok(0);
        };
        Ok(TryInto::<u8>::try_into(
            self.mapq
                .iter()
                .map(|m| u32::from(*m))
                .sum::<u32>()
                .div_ceil(length),
        )?)
    }
}

// https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009
pub fn get_aligned_pairs(
    cg: impl Iterator<Item = (Kind, usize)>,
    pos: usize,
    min_ins_size: usize,
    min_del_size: usize,
) -> eyre::Result<Vec<(usize, usize, Kind)>> {
    let mut pos: usize = pos;
    let mut qpos: usize = 0;
    let mut pairs = vec![];
    // Matches only
    for (op, l) in cg {
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                    qpos += 1
                }
                pos += l
            }
            Kind::Pad => {
                qpos += l;
                continue;
            }
            // Track indels and softclips.
            // Ignore small indels.
            Kind::Insertion | Kind::SoftClip => {
                if op == Kind::Insertion && l < min_ins_size {
                    for _ in pos..(pos + l) {
                        qpos += 1
                    }
                } else {
                    pairs.push((qpos, pos, op));
                    qpos += l
                }
            }
            Kind::Deletion => {
                if op == Kind::Deletion && l < min_del_size {
                    for i in pos..(pos + l) {
                        pairs.push((qpos, i, Kind::Match));
                        qpos += 1
                    }
                    pos += l
                } else {
                    for i in pos..(pos + l) {
                        pairs.push((qpos, i, op));
                    }
                    pos += l
                }
            }
            Kind::HardClip => {
                continue;
            }
            Kind::Skip => pos += l,
        }
    }
    Ok(pairs)
}

macro_rules! pileup {
    ($read:ident, $aln_pairs:ident, $st:ident, $end:ident, $pileup_infos:ident) => {
        // If within region of interest.
        for (_qpos, refpos, kind) in $aln_pairs
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= $st && *refpos <= $end)
        {
            let pos = refpos - $st;
            let pileup_info = &mut $pileup_infos[pos];

            match kind {
                Kind::Deletion | Kind::Insertion => {
                    pileup_info.n_indel += 1;
                    continue;
                }
                Kind::SoftClip => {
                    pileup_info.n_softclip += 1;
                    continue;
                }
                Kind::SequenceMismatch => {
                    pileup_info.n_mismatch += 1;
                }
                _ => (),
            }
            pileup_info
                .mapq
                .push($read.mapping_quality().unwrap().get());
            pileup_info.n_cov += 1;
        }
    };
}

impl AlignmentFile {
    pub fn new(aln: impl AsRef<Path> + Debug) -> eyre::Result<Self> {
        if aln
            .as_ref()
            .extension()
            .and_then(OsStr::to_str)
            .eq(&Some("cram"))
        {
            Ok(Self::Cram(
                cram::io::indexed_reader::Builder::default()
                    .build_from_path(&aln)
                    .with_context(|| format!("Cannot read indexed CRAM file ({aln:?})"))?,
            ))
        } else {
            Ok(Self::Bam(
                bam::io::indexed_reader::Builder::default()
                    .build_from_path(&aln)
                    .with_context(|| format!("Cannot read indexed BAM file ({aln:?})"))?,
            ))
        }
    }
    pub fn header(&mut self) -> eyre::Result<Header> {
        match self {
            AlignmentFile::Cram(indexed_reader) => Ok(indexed_reader.read_header()?),
            AlignmentFile::Bam(indexed_reader) => Ok(indexed_reader.read_header()?),
        }
    }

    pub fn pileup(
        &mut self,
        itv: &Interval<String>,
        min_ins_size: usize,
        min_del_size: usize,
        min_aln_length: usize,
    ) -> eyre::Result<PileupSummary> {
        let st = TryInto::<usize>::try_into(itv.first)?.clamp(1, usize::MAX);
        let end: usize = itv.last.try_into()?;
        let length = itv.len();
        // Query entire contig.
        let region = Region::new(
            &*itv.metadata,
            Position::try_from(st)?..=Position::try_from(end)?,
        );
        log::info!("Generating pileup over {}:{st}-{end}.", region.name());

        let mut pileup_infos: Vec<PileupInfo> = vec![PileupInfo::default(); length.try_into()?];
        // Reduce some redundancy with macro.
        // https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignmentfile.pyx#L1458
        match self {
            AlignmentFile::Cram(indexed_reader) => {
                let header: noodles::sam::Header = indexed_reader.read_header()?;
                let query: cram::io::reader::Query<'_, File> =
                    indexed_reader.query(&header, &region)?;
                for rec in query
                    .into_iter()
                    .flatten()
                    .filter(|aln| aln.sequence().len() > min_aln_length)
                {
                    let cg: &noodles::sam::alignment::record_buf::Cigar = rec.cigar();
                    let aln_pairs = get_aligned_pairs(
                        cg.iter().flatten().map(|op| (op.kind(), op.len())),
                        rec.alignment_start().unwrap().get(),
                        min_ins_size,
                        min_del_size,
                    )?;
                    pileup!(rec, aln_pairs, st, end, pileup_infos)
                }
            }
            AlignmentFile::Bam(indexed_reader) => {
                let header: noodles::sam::Header = indexed_reader.read_header()?;
                let query: bam::io::reader::Query<'_, bgzf::Reader<File>> =
                    indexed_reader.query(&header, &region)?;
                for rec in query
                    .into_iter()
                    .flatten()
                    .filter(|aln| aln.sequence().len() > min_aln_length)
                {
                    let cg: bam::record::Cigar<'_> = rec.cigar();
                    let aln_pairs = get_aligned_pairs(
                        cg.iter().flatten().map(|op| (op.kind(), op.len())),
                        rec.alignment_start().unwrap()?.get(),
                        min_ins_size,
                        min_del_size,
                    )?;
                    pileup!(rec, aln_pairs, st, end, pileup_infos)
                }
            }
        }
        log::info!("Finished pileup over {}:{st}-{end}.", region.name());

        Ok(PileupSummary {
            region,
            pileups: pileup_infos,
        })
    }
}

pub(crate) fn merge_pileup_info(
    pileup: Vec<PileupInfo>,
    st: u64,
    end: u64,
    cfg: &Config,
) -> eyre::Result<DataFrame> {
    let (
        mut cov_cnts,
        mut mismatch_cnts,
        mut mapq_mean_cnts,
        mut mapq_max_cnts,
        mut indel_cnts,
        mut softclip_cnts,
    ) = (
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
    );
    // Choose pileup function.
    // IMPORTANT: For false duplication detection, we need to be absolutely sure since we only have coverage and mapq to go off of.
    // * Max is generally best here as we only care if one read is high mapq.
    let pileup_fn: Box<dyn Fn(&PileupInfo) -> u8> = Box::new(match cfg.mapq.mapq_agg_fn {
        PileupMAPQFn::Mean => |p: &PileupInfo| p.mean_mapq().unwrap_or_default(),
        PileupMAPQFn::Median => |p: &PileupInfo| p.median_mapq().unwrap_or_default(),
    });
    for p in pileup.into_iter() {
        cov_cnts.push(p.n_cov);
        mismatch_cnts.push(p.n_mismatch);
        mapq_max_cnts.push(p.mapq.iter().max().cloned().unwrap_or_default());
        mapq_mean_cnts.push(pileup_fn(&p));
        indel_cnts.push(p.n_indel);
        softclip_cnts.push(p.n_softclip);
    }
    let mut lf = DataFrame::new(vec![
        Column::new("pos".into(), st..end + 1),
        Column::new("cov".into(), cov_cnts),
        Column::new("mismatch".into(), mismatch_cnts),
        Column::new("mapq_max".into(), mapq_max_cnts),
        Column::new("mapq".into(), mapq_mean_cnts),
        Column::new("indel".into(), indel_cnts),
        Column::new("softclip".into(), softclip_cnts),
    ])?
    .lazy();

    for (colname, window_size) in [
        ("cov", cfg.cov.rolling_mean_window),
        ("mismatch", cfg.mismatch.rolling_mean_window),
        ("indel", cfg.indel.rolling_mean_window),
    ] {
        if let Some(window_size) = window_size {
            lf = lf.with_column(col(colname).rolling_mean(RollingOptionsFixedWindow {
                window_size,
                center: true,
                ..Default::default()
            }))
        };
    }
    Ok(lf.collect()?)
}

#[cfg(test)]
mod test {
    use crate::{
        config::Config,
        pileup::{merge_pileup_info, AlignmentFile, PileupInfo, PileupSummary},
    };
    use noodles::core::{Position, Region};
    use polars::df;

    #[test]
    fn test_pileup() {
        let mut bam = AlignmentFile::new("test/pileup/input/test.bam").unwrap();
        let itv = coitrees::Interval::new(
            9667238,
            9667240,
            "K1463_2281_chr15_contig-0000423".to_owned(),
        );
        let res = bam.pileup(&itv, 1, 1, 0).unwrap();
        assert_eq!(
            res,
            PileupSummary {
                region: Region::new(
                    "K1463_2281_chr15_contig-0000423",
                    Position::new(9667238).unwrap()..=Position::new(9667240).unwrap()
                ),
                pileups: [
                    PileupInfo {
                        n_cov: 41,
                        n_mismatch: 0,
                        n_indel: 40,
                        n_softclip: 0,
                        mapq: [
                            60, 60, 60, 60, 60, 18, 34, 60, 35, 60, 60, 33, 30, 60, 33, 34, 33, 31,
                            33, 36, 32, 32, 60, 35, 33, 36, 31, 35, 35, 33, 33, 34, 35, 60, 33, 60,
                            60, 60, 60, 60, 60
                        ]
                        .to_vec()
                    },
                    PileupInfo {
                        n_cov: 41,
                        n_mismatch: 0,
                        n_indel: 0,
                        n_softclip: 0,
                        mapq: [
                            60, 60, 60, 60, 60, 18, 34, 60, 35, 60, 60, 33, 30, 60, 33, 34, 33, 31,
                            33, 36, 32, 32, 60, 35, 33, 36, 31, 35, 35, 33, 33, 34, 35, 60, 33, 60,
                            60, 60, 60, 60, 60
                        ]
                        .to_vec()
                    },
                    PileupInfo {
                        n_cov: 41,
                        n_mismatch: 0,
                        n_indel: 38,
                        n_softclip: 0,
                        mapq: [
                            60, 60, 60, 60, 60, 18, 34, 60, 35, 60, 60, 33, 30, 60, 33, 34, 33, 31,
                            33, 36, 32, 32, 60, 35, 33, 36, 31, 35, 35, 33, 33, 34, 35, 60, 33, 60,
                            60, 60, 60, 60, 60
                        ]
                        .to_vec()
                    }
                ]
                .to_vec()
            }
        );
    }

    #[test]
    fn test_pileup_summary_df() {
        let mut bam = AlignmentFile::new("test/pileup/input/test.bam").unwrap();
        let itv = coitrees::Interval::new(
            9667238,
            9667240,
            "K1463_2281_chr15_contig-0000423".to_owned(),
        );
        let res = bam.pileup(&itv, 1, 1, 0).unwrap();

        let config = Config::default();
        let df_pileup =
            merge_pileup_info(res.pileups, itv.first as u64, itv.last as u64, &config).unwrap();
        assert_eq!(
            df_pileup,
            df!(
                "pos" => [9667238, 9667239, 9667240],
                "cov" => [41;3],
                "mismatch" => [0; 3],
                "mapq_max" => [60; 3],
                "mapq" => [45; 3],
                "indel" => [20.0, 26.0, 19.0],
                "softclip" => [0; 3],
            )
            .unwrap()
        );
    }
}
