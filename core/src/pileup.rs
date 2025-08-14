use coitrees::{GenericInterval, Interval};
use eyre::{bail, Context};
use itertools::Itertools;
use noodles::{
    bam::{self},
    bgzf,
    core::{Position, Region},
    cram,
    sam::{
        alignment::{
            record::{
                cigar::{op::Kind, Op},
                data::field::Tag,
                Cigar, Flags,
            },
            record_buf::{data::field::Value, Data},
        },
        Header,
    },
};
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::{ffi::OsStr, fmt::Debug, fs::File, path::Path};

use crate::config::Config;

#[derive(Debug, PartialEq, Eq, Clone)]
enum DiffToken {
    Identical,
    Substitution,
    Insertion,
    Deletion,
    Intron,
    Base,
    Number,
    Invalid(u8),
}

impl From<u8> for DiffToken {
    fn from(value: u8) -> Self {
        match value {
            b'=' | b':' => DiffToken::Identical,
            b'*' => DiffToken::Substitution,
            b'+' => DiffToken::Insertion,
            b'-' => DiffToken::Deletion,
            b'~' => DiffToken::Intron,
            b'0'..=b'9' => DiffToken::Number,
            b'a' | b't' | b'g' | b'c' | b'n' => DiffToken::Base,
            _ => DiffToken::Invalid(value),
        }
    }
}

impl TryFrom<DiffToken> for Kind {
    type Error = eyre::Report;

    fn try_from(value: DiffToken) -> Result<Self, Self::Error> {
        Ok(match value {
            DiffToken::Identical => Kind::SequenceMatch,
            DiffToken::Substitution => Kind::SequenceMismatch,
            DiffToken::Insertion => Kind::Insertion,
            DiffToken::Deletion => Kind::Deletion,
            _ => bail!("Cannot convert {value:?} to Kind."),
        })
    }
}

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

fn update_cigar<'a>(
    cg: impl Iterator<Item = &'a Op>,
    tags: &Data,
) -> eyre::Result<Option<Vec<Op>>> {
    // Is not correct cigar format. Try to convert if cs provided.
    if cg.into_iter().any(|op| op.kind() == Kind::Match) {
        // :10*at:5-ac:6
        if let Some(Value::String(difference_string)) = tags.get(&Tag::new(b'c', b's')) {
            let mut curr_op: Option<DiffToken> = None;
            let mut new_ops = vec![];
            for (tk, elems) in &difference_string
                .iter()
                .chunk_by(|c| Into::<DiffToken>::into(**c))
            {
                let elems: Vec<&u8> = elems.collect();

                match (curr_op.as_ref(), &tk) {
                    (None, DiffToken::Identical)
                    | (None, DiffToken::Substitution)
                    | (None, DiffToken::Insertion)
                    | (None, DiffToken::Deletion)
                    | (None, DiffToken::Intron) => curr_op = Some(tk.clone()),
                    (None, _) => bail!("Invalid starting token: {tk:?}"),
                    (Some(op), DiffToken::Base) => {
                        let op_kind: Kind = op.clone().try_into()?;
                        new_ops.push(Op::new(op_kind, elems.len()));
                        curr_op.take();
                    }
                    (Some(op), DiffToken::Number) => {
                        let op_kind: Kind = op.clone().try_into()?;
                        new_ops.push(Op::new(
                            op_kind,
                            elems.into_iter().map(|e| *e as char).join("").parse()?,
                        ));
                        curr_op.take();
                    }
                    (Some(op), DiffToken::Identical)
                    | (Some(op), DiffToken::Substitution)
                    | (Some(op), DiffToken::Insertion)
                    | (Some(op), DiffToken::Deletion)
                    | (Some(op), DiffToken::Intron)
                    | (Some(op), DiffToken::Invalid(_)) => bail!(
                        "Invalid matching token: {op:?}: {}",
                        elems.into_iter().map(|e| *e as char).join("")
                    ),
                }
            }
            Ok(Some(new_ops))
        } else {
            bail!("Invalid cigar string. Must be extended cigar (=/X) or include cs tag.")
        }
    } else {
        Ok(None)
    }
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
                for rec in query.into_iter().flatten().filter(|aln| {
                    aln.sequence().len() > min_aln_length && !aln.flags().contains(Flags::SECONDARY)
                }) {
                    let cg: &noodles::sam::alignment::record_buf::Cigar = rec.cigar();
                    // Update cigar if cs tag is available. MD is annoying so not implemented.
                    let aln_pairs =
                        if let Some(updated_cg) = update_cigar(cg.as_ref().iter(), rec.data())? {
                            get_aligned_pairs(
                                updated_cg.into_iter().map(|op| (op.kind(), op.len())),
                                rec.alignment_start().unwrap().get(),
                                min_ins_size,
                                min_del_size,
                            )?
                        } else {
                            get_aligned_pairs(
                                cg.iter().flatten().map(|op| (op.kind(), op.len())),
                                rec.alignment_start().unwrap().get(),
                                min_ins_size,
                                min_del_size,
                            )?
                        };
                    pileup!(rec, aln_pairs, st, end, pileup_infos)
                }
            }
            AlignmentFile::Bam(indexed_reader) => {
                let header: noodles::sam::Header = indexed_reader.read_header()?;
                let query: bam::io::reader::Query<'_, bgzf::Reader<File>> =
                    indexed_reader.query(&header, &region)?;
                for rec in query.into_iter().flatten().filter(|aln| {
                    aln.sequence().len() > min_aln_length && !aln.flags().contains(Flags::SECONDARY)
                }) {
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
        pileup::{merge_pileup_info, update_cigar, AlignmentFile, PileupInfo, PileupSummary},
    };
    use noodles::core::{Position, Region};
    use noodles::sam::alignment::{
        record::{
            cigar::{op::Kind, Op},
            data::field::Tag,
        },
        record_buf::{data::field::Value, Data},
    };
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
                "pos" => [9667238u64, 9667239u64, 9667240u64],
                "cov" => [41u32;3],
                "mismatch" => [0u32; 3],
                "mapq_max" => [60u8; 3],
                "mapq" => [45u8; 3],
                "indel" => [40u32, 0u32, 38u32],
                "softclip" => [0u32; 3],
            )
            .unwrap()
        );
    }

    fn cigar_data() -> (Vec<Op>, Data, [Op; 5]) {
        let (tag, value) = (Tag::new(b'c', b's'), Value::from(":10*at:5-ac:6"));
        let data: Data = [(tag, value.clone())].into_iter().collect();

        let cg = vec![
            Op::new(Kind::Match, 16),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 6),
        ];
        const EXPECTED: [Op; 5] = [
            Op::new(Kind::SequenceMatch, 10),
            Op::new(Kind::SequenceMismatch, 2),
            Op::new(Kind::SequenceMatch, 5),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::SequenceMatch, 6),
        ];
        (cg, data, EXPECTED)
    }
    #[test]
    fn test_update_cigar() {
        let (cg, data, cg_exp) = cigar_data();
        let new_cg = update_cigar(cg.iter(), &data).unwrap();

        assert_eq!(new_cg.unwrap(), cg_exp)
    }

    #[test]
    fn test_update_cigar_no_change() {
        let (_, data, cg_exp) = cigar_data();
        // Already proper cigar
        let new_cg = update_cigar(cg_exp.iter(), &data).unwrap();
        assert!(new_cg.is_none())
    }

    #[test]
    #[should_panic]
    fn test_update_cigar_no_cs_tag() {
        let (cg, mut data, _) = cigar_data();
        data.clear();
        update_cigar(cg.iter(), &data).unwrap();
    }
}
