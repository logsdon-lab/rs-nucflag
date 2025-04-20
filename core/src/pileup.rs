use itertools::Itertools;
use noodles::{
    bam, bgzf, cram,
    core::{Region, Position},
    fasta::{self, repository},
    sam::{
        alignment::record::{cigar::op::Kind, Cigar},
        Header,
    },
};
use std::{fs::File, path::Path};
use coitrees::{GenericInterval, Interval};


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
pub fn get_aligned_pairs(
    cg: impl Iterator<Item = (Kind, usize)>,
    pos: usize,
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

macro_rules! pileup {
    ($read:ident, $aln_pairs:ident, $st:ident, $end:ident, $pileup_infos:ident) => {
        // If within region of interest.
        for (_qpos, refpos, kind) in $aln_pairs
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= $st && *refpos <= $end)
        {
            let pos = refpos - $st;
            let pileup_info = &mut $pileup_infos[pos];
            pileup_info
                .mapq
                .push($read.mapping_quality().unwrap().get());

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

            if $read.flags().is_supplementary() {
                pileup_info.n_supp += 1
            }
        }
    };
}

impl AlignmentFile {
    pub fn new(aln: impl AsRef<Path>, fasta: Option<impl AsRef<Path>>) -> eyre::Result<Self> {
        // Try to build cram.
        let is_cram = cram::io::indexed_reader::Builder::default()
            .build_from_path(aln.as_ref())
            .ok();
        // Then check if fasta and cram.
        if let Some(fasta) = is_cram.and(fasta) {
            let fasta_fh = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;
            let reference_sequence_repository =
                fasta::Repository::new(repository::adapters::IndexedReader::new(fasta_fh));
            Ok(Self::Cram(
                cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(reference_sequence_repository)
                    .build_from_path(aln)?,
            ))
        } else {
            // Otherwise, assume bam.
            Ok(Self::Bam(
                bam::io::indexed_reader::Builder::default().build_from_path(&aln)?,
            ))
        }
    }
    pub fn header(&mut self) -> eyre::Result<Header> {
        match self {
            AlignmentFile::Cram(indexed_reader) => Ok(indexed_reader.read_header()?),
            AlignmentFile::Bam(indexed_reader) => Ok(indexed_reader.read_header()?),
        }
    }

    pub fn pileup(&mut self, itv: &Interval<String>) -> eyre::Result<PileupSummary> {
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
                for read in query.into_iter().flatten() {
                    let cg: &noodles::sam::alignment::record_buf::Cigar = read.cigar();
                    let aln_pairs = get_aligned_pairs(
                        cg.iter().flatten().map(|op| (op.kind(), op.len())),
                        read.alignment_start().unwrap().get(),
                    )?;
                    pileup!(read, aln_pairs, st, end, pileup_infos)
                }
            }
            AlignmentFile::Bam(indexed_reader) => {
                let header: noodles::sam::Header = indexed_reader.read_header()?;
                let query: bam::io::reader::Query<'_, bgzf::Reader<File>> =
                    indexed_reader.query(&header, &region)?;
                for read in query.into_iter().flatten() {
                    let cg: bam::record::Cigar<'_> = read.cigar();
                    let aln_pairs = get_aligned_pairs(
                        cg.iter().flatten().map(|op| (op.kind(), op.len())),
                        read.alignment_start().unwrap()?.get(),
                    )?;
                    pileup!(read, aln_pairs, st, end, pileup_infos)
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

#[cfg(test)]
mod test {
    use crate::pileup::{AlignmentFile, PileupInfo, PileupSummary};
    use noodles::core::{Position, Region};

    #[test]
    fn test_pileup() {
        let mut bam = AlignmentFile::new("test/pileup/test.bam", None::<&str>).unwrap();
        let itv = coitrees::Interval::new(
            9667238,
            9667240,
            "K1463_2281_chr15_contig-0000423".to_owned(),
        );
        let res = bam.pileup(&itv).unwrap();
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
