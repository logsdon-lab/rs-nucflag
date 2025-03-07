use std::fs::File;

use coitrees::{GenericInterval, Interval};
use eyre::ContextCompat;
use itertools::Itertools;
use noodles::{
    bam,
    core::{Position, Region},
    sam::alignment::record::cigar::op::Kind,
};

#[derive(Debug, Clone, Default)]
pub struct PileupInfo {
    pub n_a: u64,
    pub n_t: u64,
    pub n_g: u64,
    pub n_c: u64,
    pub n_sec: u64,
    pub n_sup: u64,
    pub mapq: Vec<u8>,
}

pub struct PileupSummary {
    pub region: Region,
    pub pileups: Vec<PileupInfo>,
}

pub struct DepthSummary {
    pub region: Region,
    pub avg_depth: u64,
    // pub depth: Vec<u16>,
}

impl PileupInfo {
    pub fn median_mapq(&self) -> Option<u8> {
        self.mapq.iter().sorted().nth(self.mapq.len() / 2).cloned()
    }
    #[allow(unused)]
    pub fn mean_mapq(&self) -> eyre::Result<u8> {
        self.mapq
            .iter()
            .sum::<u8>()
            .checked_sub(TryInto::<u8>::try_into(self.mapq.len())?)
            .context("No mapq.")
    }
    pub fn counts(&self) -> impl Iterator<Item = u64> {
        [self.n_a, self.n_t, self.n_g, self.n_c]
            .into_iter()
            .sorted_by(|a, b| a.cmp(b))
            .rev()
    }
    pub fn first(&self) -> u64 {
        // Will never be empty so unwrap safe.
        self.counts().next().unwrap()
    }

    pub fn second(&self) -> u64 {
        // Will never be empty so unwrap safe.
        self.counts().nth(1).unwrap()
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
            Kind::Insertion | Kind::SoftClip | Kind::Pad => qpos += l,
            Kind::Deletion => pos += l,
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
    let st: usize = itv.first.try_into()?;
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
        let seq = read.sequence();
        let mapq = read.mapping_quality().unwrap().get();
        let flags = read.flags();
        // If within region of interest.
        for (qpos, refpos, _) in get_aligned_pairs(&read)?
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
        {
            let pos = refpos - st;
            let pileup_info = &mut pileup_infos[pos];
            pileup_info.mapq.push(mapq);

            if flags.is_supplementary() {
                pileup_info.n_sup += 1
            }
            if flags.is_secondary() {
                pileup_info.n_sec += 1
            }
            let bp = seq.get(qpos).unwrap();
            match bp {
                b'A' => pileup_info.n_a += 1,
                b'C' => pileup_info.n_c += 1,
                b'G' => pileup_info.n_g += 1,
                b'T' => pileup_info.n_t += 1,
                b'N' => (),
                _ => log::debug!("Character not recognized at pos {pos}: {bp:?}"),
            }
        }
    }
    log::info!("Finished pileup over {}:{st}-{end}.", region.name());

    Ok(PileupSummary {
        region,
        pileups: pileup_infos,
    })
}
