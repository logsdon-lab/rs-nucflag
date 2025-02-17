use std::{
    fs::File,
    io::{BufWriter, Write},
    str::FromStr,
};

use coitrees::{self as ct, GenericInterval, Interval};
use draw::draw_nucfreq;
use eyre::bail;
use io::write_tsv;
use itertools::Itertools;
use noodles::{
    bam::{self},
    core::{Position, Region},
    sam::alignment::record::cigar::op::Kind,
};
use peak::{find_peaks, merge_peaks, Peak};
use polars::prelude::*;

mod draw;
mod intervals;
mod io;
mod misassembly;
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

fn peak_cols<'a>(
    pos: &'a Column,
    peak_col: &'a Column,
) -> eyre::Result<impl Iterator<Item = (u64, Peak)> + use<'a>> {
    Ok(pos
        .u64()?
        .iter()
        .flatten()
        .zip(
            peak_col
                .str()?
                .iter()
                .flat_map(|p| Peak::from_str(p.unwrap_or_default())),
        )
        .flat_map(|(pos, peak)| (peak != Peak::Null).then_some((pos, peak))))
}

fn df_pileup_info(pileup: Vec<PileupInfo>, st: u64, end: u64) -> eyre::Result<DataFrame> {
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
        mapq.push(p.median_mapq().unwrap_or(0));
        sec.push(p.n_sec);
        sup.push(p.n_sup);
    }
    DataFrame::new(vec![
        Column::new("pos".into(), st..end + 1),
        Column::new("first".into(), first),
        Column::new("second".into(), second),
        Column::new("mapq".into(), mapq),
        Column::new("sec".into(), sec),
        Column::new("sup".into(), sup),
    ])
    .map_err(Into::into)
}

// https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
fn main() -> eyre::Result<()> {
    let window_size = 500_000;
    let n_stdevs = 5;
    let thr_second_perc = 0.25;
    let dim = (4000, 700);
    let output_png = "out.png";

    let (st, end) = (3021508u64, 8691473u64);
    let ctg = "haplotype2-0000133".to_owned();

    let file = "test/standard/aln.bam";
    let mut bam = bam::io::indexed_reader::Builder::default().build_from_path(file)?;
    let pileup = pileup(
        &mut bam,
        &ct::Interval::new(st.try_into()?, end.try_into()?, ctg.to_owned()),
    )?;

    let df_pileup_info = df_pileup_info(pileup, st.try_into()?, end.try_into()?)?;
    let df_original_data = df_pileup_info.select(["first", "second", "mapq"])?;
    let [first, second, mapq] = df_original_data.get_columns() else {
        bail!("Invalid number of columns. Developer error.")
    };

    let mut df = find_peaks(df_pileup_info, window_size, n_stdevs, thr_second_perc)?;
    // obed
    write_tsv(&mut df, "out.tsv")?;

    let [pos, first_peak, second_peak] = df.get_columns() else {
        bail!("Invalid number of columns. Developer error.")
    };

    let first_peaks: Vec<Interval<Peak>> =
        merge_peaks(peak_cols(pos, first_peak)?, Some(5000), Some(100))?;

    let second_peaks: Vec<Interval<Peak>> =
        merge_peaks(peak_cols(pos, second_peak)?, Some(5000), Some(100))?;

    let bed_first = File::create("first.bed")?;
    let mut bed_first_writer = BufWriter::new(bed_first);
    for peak in first_peaks.iter() {
        writeln!(
            &mut bed_first_writer,
            "{ctg}\t{}\t{}\t{:?}",
            peak.first, peak.last, peak.metadata
        )?;
    }

    let bed_second = File::create("second.bed")?;
    let mut bed_second_writer = BufWriter::new(bed_second);
    for peak in second_peaks.iter() {
        writeln!(
            &mut bed_second_writer,
            "{ctg}\t{}\t{}\t{:?}",
            peak.first, peak.last, peak.metadata
        )?;
    }

    draw_nucfreq(
        output_png,
        dim,
        &ctg,
        pos,
        first,
        second,
        mapq,
        first_peaks,
        second_peaks,
    )?;

    Ok(())
}
