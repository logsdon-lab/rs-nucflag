use std::fs::File;

use coitrees::{GenericInterval, Interval};
use itertools::Itertools;
use noodles::{
    bam,
    core::{Position, Region},
    sam::alignment::record::cigar::op::Kind,
};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PileupInfo {
    pub n_cov: u64,
    pub n_mismatch: u64,
    pub n_indel: u64,
    pub n_supp: u64,
    pub n_softclip: u64,
    pub mapq: Vec<u8>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct PileupSummary {
    pub region: Region,
    pub pileups: Vec<PileupInfo>,
}

impl PileupInfo {
    pub fn median_mapq(&self) -> Option<u8> {
        let length = self.mapq.len();
        let midpt = length / 2;
        if length % 2 == 0 {
            Some(
                self.mapq
                    .iter()
                    .sorted()
                    .get(midpt..=midpt + 1)
                    .sum::<u8>()
                    .div_ceil(2),
            )
        } else {
            self.mapq.iter().sorted().nth(self.mapq.len() / 2).cloned()
        }
    }
    pub fn mean_mapq(&self) -> eyre::Result<u8> {
        let Some(length) = TryInto::<u8>::try_into(self.mapq.len())
            .ok()
            .filter(|l| *l > 0)
        else {
            return Ok(0);
        };
        Ok(self.mapq.iter().sum::<u8>().div_ceil(length))
    }
}

// https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009
pub fn get_aligned_pairs(read: &bam::Record) -> eyre::Result<Vec<(usize, usize, Kind)>> {
    let cg = read.cigar();
    let mut pos: usize = read.alignment_start().unwrap()?.get();
    let mut qpos: usize = 0;
    let mut pairs = vec![];
    // Matches only
    for (op, l) in cg.iter().flatten().map(|op| (op.kind(), op.len())) {
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
            Kind::Insertion | Kind::SoftClip => {
                pairs.push((qpos, pos, op));
                qpos += l
            }
            Kind::Deletion => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                }
                pos += l
            }
            Kind::HardClip => {
                continue;
            }
            Kind::Skip => pos += l,
        }
    }
    Ok(pairs)
}

pub fn pileup(
    bam: &mut bam::io::IndexedReader<noodles::bgzf::Reader<File>>,
    itv: &Interval<String>,
) -> eyre::Result<PileupSummary> {
    let st = TryInto::<usize>::try_into(itv.first)?.clamp(1, usize::MAX);
    let end: usize = itv.last.try_into()?;
    let length: usize = itv.len().try_into()?;

    // Query entire contig.
    let region = Region::new(
        &*itv.metadata,
        Position::try_from(st)?..=Position::try_from(end)?,
    );

    let header = bam.read_header()?;

    // https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignmentfile.pyx#L1458
    let query = bam.query(&header, &region)?;
    let mut pileup_infos: Vec<PileupInfo> = vec![PileupInfo::default(); length];

    log::info!("Generating pileup over {}:{st}-{end}.", region.name());
    for read in query.into_iter().flatten() {
        let mapq = read.mapping_quality().unwrap().get();
        // If within region of interest.
        for (_qpos, refpos, kind) in get_aligned_pairs(&read)?
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
        {
            let pos = refpos - st;
            let pileup_info = &mut pileup_infos[pos];
            pileup_info.mapq.push(mapq);

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

            pileup_info.n_cov += 1;

            if read.flags().is_supplementary() {
                pileup_info.n_supp += 1
            }
        }
    }
    log::info!("Finished pileup over {}:{st}-{end}.", region.name());

    Ok(PileupSummary {
        region,
        pileups: pileup_infos,
    })
}

#[cfg(test)]
mod test {
    use crate::pileup::{PileupInfo, PileupSummary};
    use noodles::{
        bam,
        core::{Position, Region},
    };

    use super::pileup;

    #[test]
    fn test_pileup() {
        let mut bam = bam::io::indexed_reader::Builder::default()
            .build_from_path("test/pileup/test.bam")
            .unwrap();
        let itv = coitrees::Interval::new(
            9667238,
            9667240,
            "K1463_2281_chr15_contig-0000423".to_owned(),
        );
        let res = pileup(&mut bam, &itv).unwrap();
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
                        n_supp: 0,
                        n_softclip: 0,
                        mapq: [
                            60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 18, 18, 34, 34, 60, 60, 35, 35,
                            60, 60, 60, 60, 33, 33, 30, 30, 60, 60, 33, 33, 34, 34, 33, 33, 31, 31,
                            33, 33, 36, 36, 32, 32, 32, 32, 60, 60, 35, 35, 33, 33, 36, 36, 31, 35,
                            35, 35, 35, 33, 33, 33, 33, 34, 34, 35, 35, 60, 60, 33, 33, 60, 60, 60,
                            60, 60, 60, 60, 60, 60, 60, 60, 60
                        ]
                        .to_vec()
                    },
                    PileupInfo {
                        n_cov: 41,
                        n_mismatch: 0,
                        n_indel: 0,
                        n_supp: 0,
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
                        n_supp: 0,
                        n_softclip: 0,
                        mapq: [
                            60, 60, 60, 60, 60, 60, 60, 60, 18, 18, 34, 34, 60, 60, 35, 35, 60, 60,
                            60, 60, 33, 33, 30, 30, 60, 60, 33, 33, 34, 34, 33, 33, 31, 31, 33, 33,
                            36, 36, 32, 32, 32, 32, 60, 60, 35, 35, 33, 33, 36, 36, 31, 31, 35, 35,
                            35, 35, 33, 33, 33, 33, 34, 34, 35, 35, 60, 60, 33, 33, 60, 60, 60, 60,
                            60, 60, 60, 60, 60, 60, 60
                        ]
                        .to_vec()
                    }
                ]
                .to_vec()
            }
        );
    }
}
