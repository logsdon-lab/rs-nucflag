use std::{path::Path, str::FromStr};

use crate::{
    config::Config,
    draw::draw_nucfreq,
    io::write_tsv,
    misassembly::MisassemblyType,
    peak::{find_peaks, merge_positions},
    pileup::{pileup, PileupInfo},
    Interval,
};
use itertools::multizip;
use noodles::bam::{self};
use polars::prelude::*;

fn merge_pileup_info(pileup: Vec<PileupInfo>, st: u64, end: u64) -> eyre::Result<DataFrame> {
    let (mut first_cnts, mut second_cnts, mut mapq_cnts, mut sec_cnts, mut sup_cnts) = (
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
    );
    for p in pileup.into_iter() {
        first_cnts.push(p.first());
        second_cnts.push(p.second());
        mapq_cnts.push(p.median_mapq().unwrap_or(0));
        sec_cnts.push(p.n_sec);
        sup_cnts.push(p.n_sup);
    }
    DataFrame::new(vec![
        Column::new("pos".into(), st..end + 1),
        Column::new("first".into(), first_cnts),
        Column::new("second".into(), second_cnts),
        Column::new("mapq".into(), mapq_cnts),
        Column::new("sec".into(), sec_cnts),
        Column::new("sup".into(), sup_cnts),
    ])
    .map_err(Into::into)
}

fn build_interval(
    positions: &[(u64, (MisassemblyType, u64))],
) -> eyre::Result<Interval<(MisassemblyType, u64)>> {
    let (st, end, (typ, cov)) =
        positions
            .iter()
            .fold((u64::MAX, 0u64, (MisassemblyType::Null, 0)), |acc, x| {
                let (st, end, (_, cov_sum)) = acc;
                let (curr_pos, (curr_typ, curr_cov)) = x;
                (
                    std::cmp::min(st, *curr_pos),
                    std::cmp::max(end, *curr_pos),
                    (curr_typ.clone(), curr_cov + cov_sum),
                )
            });
    Ok(Interval::new(st, end, (typ, cov)))
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
) -> eyre::Result<(String, Vec<Interval<(MisassemblyType, u64)>>)> {
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.st, itv.end);
    let mut bam = bam::io::indexed_reader::Builder::default().build_from_path(&bamfile)?;
    let pileup = pileup(&mut bam, itv)?;

    let df_raw_pileup = merge_pileup_info(pileup, st, end)?;

    // TODO: Figure out both collapses, misjoins, and false duplications
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");

    // This is never filtered so safe to left join without removing data.
    let lf_first_peaks = find_peaks(
        df_raw_pileup.select(["pos", "first"])?,
        cfg.first.window_size,
        cfg.first.n_zscores_collapse,
        None,
    )?;
    let lf_second_peaks = find_peaks(
        df_raw_pileup.select(["pos", "second"])?,
        cfg.second.window_size,
        cfg.second.n_zscores,
        Some(cfg.second.min_perc),
    )?;

    let mut df_pileup = lf_first_peaks
        .join(
            lf_second_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .join(
            df_raw_pileup.select(["pos", "mapq"])?.lazy(),
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        // Fill nulls from join.
        .with_column(col("second").fill_null(lit(0)))
        .with_columns([
            col("second_peak").fill_null(lit("null")),
            // Calculate het ratio.
            (col("second").cast(DataType::Float32)
                / (col("first") + col("second")).cast(DataType::Float32))
            .alias("het_ratio"),
        ])
        .with_columns([
            // collapse_other
            // Regions with higher than expected het ratio.
            when(
                col("het_ratio")
                    .gt(lit(cfg.second.thr_het_ratio))
                    .and(col("second_peak").eq(lit("high")))
                    .and(col("first_peak").eq(lit("null"))),
            )
            .then(lit("collapse_other"))
            // misjoin
            // Regions with either:
            // * Zero coverage. Might be scaffold or misjoined contig. Without reference, not known.
            // * A dip in coverage with a high het ratio.
            .when(
                col("first").eq(lit(0)).or(col("first_peak")
                    .eq(lit("low"))
                    .and(col("second_peak").eq(lit("high")))),
            )
            .then(lit("misjoin"))
            // collapse_var
            // Regions with high coverage and a high coverage of another variant.
            .when(
                col("second_peak")
                    .eq(lit("high"))
                    .and(col("first_peak").eq(lit("high"))),
            )
            .then(lit("collapse_var"))
            // collapse
            // Perfect collapse of repetitive region with few to no secondary variants.
            .when(
                col("first_peak")
                    .eq(lit("high"))
                    .and(col("second_peak").eq(lit("null"))),
            )
            .then(lit("collapse"))
            // false_dupe
            // Region with half of the expected coverage and zero-mapq due to multi-mapping.
            // Either a duplicated contig or duplicated region.
            .when(
                col("first")
                    .lt_eq(col("first").median() / lit(2))
                    .and(col("mapq").eq(lit(0))),
            )
            .then(lit("false_dupe"))
            .otherwise(lit("good"))
            .alias("status"),
        ])
        .select([
            col("pos"),
            col("first"),
            col("second"),
            col("first_mean"),
            col("second_mean"),
            col("mapq"),
            col("status"),
        ])
        .collect()?;

    // Output raw coverage.
    if let Some(cov_bed) = cov_bed {
        write_tsv(&mut df_pileup, Some(cov_bed))?;
    }

    // Construct intervals.
    // Store [st,end,type,cov]
    // TODO: Must separate out for first and second.
    let itvs_all: Vec<Interval<(MisassemblyType, u64)>> = merge_positions(
        multizip((
            df_pileup.column("pos")?.u64()?.iter().flatten(),
            df_pileup.column("first")?.u64()?.iter().flatten(),
            df_pileup.column("second")?.u64()?.iter().flatten(),
            df_pileup
                .column("status")?
                .str()?
                .iter()
                .flatten()
                .map(String::from),
        ))
        .map(|(pos, first, second, status)| {
            (
                pos,
                (MisassemblyType::from_str(&status).unwrap(), first + second),
            )
        }),
        |(_pos1, (typ1, _cov1)), (_pos2, (typ2, _cov2))| typ1 == typ2,
        build_interval,
        Some(cfg.second.bp_merge.try_into()?),
        Some(cfg.second.min_bp_peak.try_into()?),
    )?;

    let misasemblies: Vec<Interval<&MisassemblyType>> = itvs_all
        .iter()
        .flat_map(|itv| {
            (itv.metadata.0 != MisassemblyType::Null).then_some(Interval::new(
                itv.st,
                itv.end,
                &itv.metadata.0,
            ))
        })
        .collect();

    log::info!(
        "Detected {} misassemblies for {ctg}:{st}-{end}.",
        misasemblies.len()
    );
    if let Some(plot_fname) = plot {
        log::info!("Plotting misassemblies for {ctg}:{st}-{end}.");
        draw_nucfreq(
            plot_fname,
            cfg.plot.dim,
            &ctg,
            &df_pileup,
            misasemblies.into_iter(),
        )?;
    }
    // Add ctg.
    Ok((ctg, itvs_all))
}
