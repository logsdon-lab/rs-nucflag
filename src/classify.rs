use std::{path::Path, str::FromStr};

use crate::{
    config::Config,
    draw::draw_nucfreq,
    io::write_tsv,
    peak::{find_peaks, merge_peaks, Peak},
    pileup::{pileup, PileupInfo},
};
use coitrees::{COITree, Interval, IntervalTree};
use eyre::bail;
use noodles::bam::{self};
use polars::prelude::*;

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

/// Classify misasemblies from alignment read coverage.
///
/// # Arguments
/// * `bamfile`: Input BAM file path. Should be indexed and filtered to only primary alignments (?).  
/// * `itv`: Interval to check.
/// * `output_bed`: Output BED file with misassemblies.
/// * `cov_dir`: Output directory to write .
/// * `plot_dir`: Output BED file with misassemblies.
pub fn classify_misassemblies(
    bamfile: impl AsRef<Path>,
    itv: Interval<String>,
    cov_bed: Option<impl AsRef<Path>>,
    plot: Option<impl AsRef<Path>>,
    cfg: Config,
) -> eyre::Result<LazyFrame> {
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.first.try_into()?, itv.last.try_into()?);
    let mut bam = bam::io::indexed_reader::Builder::default().build_from_path(&bamfile)?;
    let pileup = pileup(&mut bam, itv)?;

    let df_pileup_info = df_pileup_info(pileup, st, end)?;
    let df_original_data = df_pileup_info.select(["pos", "first", "second", "mapq"])?;

    // TODO: Figure out both collapses, misjoins, and false duplications
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");
    let lf_first_peaks = find_peaks(
        df_pileup_info.select(["pos", "first"])?,
        cfg.first.window_size,
        cfg.first.n_zscores_collapse,
        None,
    )?;
    let lf_second_peaks = find_peaks(
        df_pileup_info.select(["pos", "second"])?,
        cfg.second.window_size,
        cfg.second.n_zscores,
        Some(cfg.second.min_perc),
    )?;
    let mut df_dupes = df_pileup_info
        .select(["pos", "first", "mapq"])?
        .lazy()
        .with_column(
            (col("first")
                .lt_eq(col("first").median() / lit(2))
                .and(col("mapq").eq(0)))
            .alias("dupe"),
        )
        .collect()?;
    write_tsv(&mut df_dupes, Some("dupe.tsv"))?;

    let mut df_pileup = lf_first_peaks
        .join(
            lf_second_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .select([col("pos"), col("first_peak"), col("second_peak")])
        .with_columns([col("second_peak").fill_null(lit("none"))])
        .collect()?;

    // Output coverage.
    if let Some(cov_bed) = cov_bed {
        write_tsv(&mut df_pileup, Some(cov_bed))?;
    }

    let [pos, first_peak, second_peak] = df_pileup.get_columns() else {
        bail!("Invalid number of columns. Developer error.")
    };

    let first_peaks: Vec<Interval<Peak>> = merge_peaks(
        peak_cols(pos, first_peak)?,
        Some(cfg.first.bp_merge.try_into()?),
        Some(cfg.first.min_bp_peak.try_into()?),
    )?;

    let second_peaks: Vec<Interval<Peak>> = merge_peaks(
        peak_cols(pos, second_peak)?,
        Some(cfg.second.bp_merge.try_into()?),
        Some(cfg.second.min_bp_peak.try_into()?),
    )?;
    let first_itvs: COITree<Peak, usize> = COITree::new(&first_peaks);
    let second_itvs: COITree<Peak, usize> = COITree::new(&second_peaks);

    // TODO: Determine hets/small collapses.
    // TODO: Determine classifications.

    // let bed_first = File::create("first.bed")?;
    // let mut bed_first_writer = BufWriter::new(bed_first);
    // for peak in first_peaks.iter() {
    //     writeln!(
    //         &mut bed_first_writer,
    //         "{ctg}\t{}\t{}\t{:?}",
    //         peak.first, peak.last, peak.metadata
    //     )?;
    // }

    // let bed_second = File::create("second.bed")?;
    // let mut bed_second_writer = BufWriter::new(bed_second);
    // for peak in second_peaks.iter() {
    //     writeln!(
    //         &mut bed_second_writer,
    //         "{ctg}\t{}\t{}\t{:?}",
    //         peak.first, peak.last, peak.metadata
    //     )?;
    // }
    log::info!("Found x misassemblies for {ctg}:{st}-{end}.");
    if let Some(plot_fname) = plot {
        log::info!("Plotting misassemblies for {ctg}:{st}-{end}.");
        draw_nucfreq(
            plot_fname,
            cfg.plot.dim,
            &ctg,
            df_original_data,
            first_peaks,
            second_peaks,
        )?;
    }
    Ok(DataFrame::empty().lazy())
}
