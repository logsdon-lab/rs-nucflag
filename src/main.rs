use std::{
    fs::File,
    io::{BufWriter, Write},
    ops::Range,
};

use coitrees::{self as ct, GenericInterval, Interval};
use itertools::Itertools;
use noodles::{
    bam::{self},
    core::{Position, Region},
    sam::alignment::record::cigar::op::Kind,
};
use peak::{merge_peaks, Peak, PeaksDetector, PeaksFilter};
use plotters::prelude::*;

// https://github.com/swizard0/smoothed_z_score/blob/master/src/lib.rs
mod peak;

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

impl PileupInfo {
    fn median_mapq(&self) -> Option<u8> {
        self.mapq.iter().sorted().nth(self.mapq.len() / 2).cloned()
    }

    fn first(&self) -> u64 {
        std::cmp::max(
            std::cmp::max(self.n_a, self.n_t),
            std::cmp::max(self.n_g, self.n_c),
        )
    }

    fn second(&self) -> u64 {
        std::cmp::min(
            std::cmp::max(self.n_a, self.n_t),
            std::cmp::max(self.n_g, self.n_c),
        )
    }
}

// https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009
fn get_aligned_pairs(read: &bam::Record) -> eyre::Result<Vec<(usize, usize, Kind)>> {
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
) -> eyre::Result<Vec<PileupInfo>> {
    let st: usize = itv.first.try_into()?;
    let end: usize = itv.last.try_into()?;
    let length: usize = itv.len().try_into()?;

    let header = bam.read_header()?;
    let region = Region::new(
        &*itv.metadata,
        Position::try_from(st)?..=Position::try_from(end)?,
    );

    // https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignmentfile.pyx#L1458
    let query = bam.query(&header, &region)?;
    let mut pileup_infos: Vec<PileupInfo> = vec![PileupInfo::default(); length];

    for read in query.into_iter().flatten() {
        let seq = read.sequence();
        let mapq = read.mapping_quality().unwrap().get();
        let flags = read.flags();

        for (qpos, refpos, _) in get_aligned_pairs(&read)?
            .iter()
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
            let bp = seq.get(*qpos).unwrap();
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

    Ok(pileup_infos)
}

fn draw_nucfreq(
    outfile: &str,
    dim: (u32, u32),
    interval: &Interval<String>,
    first_data: Vec<u64>,
    second_data: Vec<u64>,
    first_peaks: Vec<Interval<Peak>>,
    second_peaks: Vec<Interval<Peak>>,
) -> eyre::Result<()> {
    let root_area = BitMapBackend::new(outfile, dim).into_drawing_area();

    root_area.fill(&WHITE)?;

    let root_area = root_area.titled(&interval.metadata, ("sans-serif", 34))?;
    let x_range: Range<u64> = interval.first.try_into()?..interval.last.try_into()?;
    let y_range = 0u64..50u64;

    let mut cc = ChartBuilder::on(&root_area)
        .margin(5)
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(x_range.clone(), y_range.clone())?;

    cc.configure_mesh()
        .x_desc("Position")
        .y_desc("Coverage")
        .draw()?;

    cc.draw_series(LineSeries::new(
        first_data
            .iter()
            .zip(x_range.clone().into_iter())
            .map(|(y, x)| (x, *y)),
        &BLACK,
    ))?;
    cc.draw_series(LineSeries::new(
        second_data
            .iter()
            .zip(x_range.clone().into_iter())
            .map(|(y, x)| (x, *y)),
        &RED,
    ))?;

    cc.draw_series(first_peaks.iter().map(|pk| {
        Rectangle::new(
            [
                (pk.first.try_into().unwrap(), y_range.end),
                (pk.last.try_into().unwrap(), y_range.start),
            ],
            ShapeStyle {
                color: BLACK.to_rgba().mix(0.5),
                filled: true,
                stroke_width: 1,
            },
        )
    }))?;

    cc.draw_series(second_peaks.iter().map(|pk| {
        Rectangle::new(
            [
                (pk.first.try_into().unwrap(), y_range.end),
                (pk.last.try_into().unwrap(), y_range.start),
            ],
            ShapeStyle {
                color: RED.to_rgba().mix(0.5),
                filled: true,
                stroke_width: 1,
            },
        )
    }))?;

    root_area.present()?;
    Ok(())
}

// https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
fn main() -> eyre::Result<()> {
    let file = "test/standard/aln.bam";
    let mut bam = bam::io::indexed_reader::Builder::default().build_from_path(file)?;

    let (st, end) = (3021508u64, 8691473u64);
    let ctg = "haplotype2-0000133".to_owned();
    let itv = ct::Interval::new(st.try_into()?, end.try_into()?, ctg.to_owned());
    let pileup = pileup(&mut bam, &itv)?;
    let (mut first, mut second, mut mapq, mut sec, mut sup) = (
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
    );
    for p in pileup.into_iter() {
        first.push(p.first());
        second.push(p.second());
        mapq.push(p.median_mapq());
        sec.push(p.n_sec);
        sup.push(p.n_sup);
    }


    // Peak calling too slow. Rewrite using polars lazy.
    let pos: Range<usize> = st.try_into()?..end.try_into()?;
    let first_peaks: Vec<Interval<Peak>> = merge_peaks(
        first
            .iter()
            .zip(pos.clone())
            .map(|pk| (pk.1, pk.0))
            .peaks(PeaksDetector::new(1_000, 10.0, 0.0), |e| *e.1 as f64)
            .map(|((i, _), p)| (i, p)),
        Some(5000),
        Some(100),
    )?;

    let second_peaks: Vec<Interval<Peak>> = merge_peaks(
        second
            .iter()
            .zip(pos)
            .filter_map(|(cov, pos)| (*cov < 2).then_some((pos, cov)))
            .peaks(PeaksDetector::new(1_000, 10.0, 0.0), |e| *e.1 as f64)
            .map(|((i, _), p)| (i, p)),
        Some(5000),
        Some(100),
    )?;

    let bed_first = File::open("first.bed")?;
    let mut bed_first_writer = BufWriter::new(bed_first);
    for peak in first_peaks.iter() {
        writeln!(
            &mut bed_first_writer,
            "{ctg}\t{}\t{}\t{:?}",
            peak.first, peak.last, peak.metadata
        )?;
    }

    let bed_second = File::open("second.bed")?;
    let mut bed_second_writer = BufWriter::new(bed_second);
    for peak in second_peaks.iter() {
        writeln!(
            &mut bed_second_writer,
            "{ctg}\t{}\t{}\t{:?}",
            peak.first, peak.last, peak.metadata
        )?;
    }

    draw_nucfreq(
        "out.png",
        (2000, 350),
        &itv,
        first,
        second,
        first_peaks,
        second_peaks,
    )?;

    Ok(())
}
