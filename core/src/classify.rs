use std::{path::Path, str::FromStr};

use crate::{
    config::Config,
    intervals::merge_overlapping_intervals,
    misassembly::MisassemblyType,
    peak::find_peaks,
    pileup::{pileup, PileupInfo},
};
use coitrees::{GenericInterval, Interval};
use itertools::{multizip, Itertools};
use noodles::bam::{self};
use polars::prelude::*;

fn merge_pileup_info(pileup: Vec<PileupInfo>, st: u64, end: u64) -> eyre::Result<DataFrame> {
    let (mut first_cnts, mut second_cnts, mut mapq_cnts) = (
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
        Vec::with_capacity(pileup.len()),
    );
    for p in pileup.into_iter() {
        first_cnts.push(p.first());
        second_cnts.push(p.second());
        mapq_cnts.push(p.median_mapq().unwrap_or(0));
    }
    DataFrame::new(vec![
        Column::new("pos".into(), st..end + 1),
        Column::new("first".into(), first_cnts),
        Column::new("second".into(), second_cnts),
        Column::new("mapq".into(), mapq_cnts),
    ])
    .map_err(Into::into)
}

fn merge_misassemblies(
    df_itvs: DataFrame,
    bp_merge: i32,
    bp_filter: i32,
    merge_across_type: bool,
) -> eyre::Result<LazyFrame> {
    let mut lf_itvs_all = df_itvs.clone().lazy();
    let df_misasm_itvs = df_itvs
        .lazy()
        .filter(col("status").neq(lit("good")))
        .collect()?;

    // Group by interval type.
    let merged_misasm_itvs: Vec<Interval<(&str, u8)>> = multizip((
        df_misasm_itvs.column("st")?.u64()?.iter().flatten(),
        df_misasm_itvs.column("end")?.u64()?.iter().flatten(),
        df_misasm_itvs.column("cov")?.u8()?.iter().flatten(),
        df_misasm_itvs.column("status")?.str()?.iter().flatten(),
    ))
    .chunk_by(|a| a.3)
    .into_iter()
    .fold(Vec::default(), |mut acc, (grp, grps)| {
        // Add base bp merge.
        acc.extend(merge_overlapping_intervals(
            grps.into_iter().map(|(st, end, cov, status)| {
                Interval::new(st as i32 - bp_merge, end as i32 + bp_merge, (status, cov))
            }),
            // Average coverage.
            |itv_1, itv_2| (itv_1.metadata.0, (itv_1.metadata.1 + itv_2.metadata.1) / 2),
            |itv: Interval<(&str, u8)>| {
                Interval::new(
                    itv.first + bp_merge,
                    itv.last - bp_merge,
                    (grp, itv.metadata.1),
                )
            },
        ));
        acc
    });

    // Merge overlapping intervals OVER status type choosing largest misassembly type.
    let final_misasm_itvs = if merge_across_type {
        merge_overlapping_intervals(
            merged_misasm_itvs.into_iter(),
            |itv_1, itv_2| {
                let largest_itv =
                    std::cmp::max_by(itv_1, itv_2, |itv_1, itv_2| itv_1.len().cmp(&itv_2.len()));
                (
                    largest_itv.metadata.0,
                    (itv_1.metadata.1 + itv_2.metadata.1) / 2,
                )
            },
            |itv| itv,
        )
    } else {
        merged_misasm_itvs
    };

    // Replace intervals between good.
    for itv in final_misasm_itvs {
        lf_itvs_all = lf_itvs_all.with_column(
            when(
                col("st")
                    .gt_eq(lit(itv.first))
                    .and(col("end").lt_eq(itv.last)),
            )
            .then(lit(itv.metadata.0))
            .otherwise(col("status"))
            .alias("status"),
        );
    }

    Ok(lf_itvs_all
        // First reduce interval groups to min/max.
        .with_column(col("status").rle_id().alias("group"))
        .group_by(["group"])
        .agg([
            col("st").min(),
            col("end").max(),
            col("cov").median(),
            col("status").first(),
        ])
        // Reset intervals smaller than filter
        .with_column(
            when((col("end") - col("st")).gt_eq(lit(bp_filter)))
                .then(col("status"))
                .otherwise(lit("good"))
                .alias("status"),
        )
        // Then reduce final interval groups to min/max.
        .with_column(col("status").rle_id().alias("group"))
        .group_by(["group"])
        .agg([
            col("st").min(),
            col("end").max(),
            col("cov").median(),
            col("status").first(),
        ])
        .sort(["st"], Default::default())
        .select([col("st"), col("end"), col("status"), col("cov")]))
}

#[derive(Debug)]
pub struct NucFlagResult {
    /// All called regions.
    pub regions: DataFrame,
    /// Pileup of regions.
    pub cov: Option<DataFrame>,
}

fn classify_peaks(
    lf_pileup: LazyFrame,
    ctg: &str,
    cfg: &Config,
    median_cov: u64
) -> eyre::Result<(DataFrame, Option<DataFrame>)> {
    let df_pileup = lf_pileup
        .with_column(
            // collapse_other
            // Regions with higher than expected het ratio.
            when(col("het_ratio").gt_eq(lit(cfg.second.thr_het_ratio)))
                .then(lit("collapse_other"))
                .otherwise(lit("good"))
                .alias("status"),
        )
        .with_column(
            // collapse_var
            // Regions with high coverage and a high coverage of another variant.
            when(
                col("first_peak")
                    .eq(lit("high"))
                    .and(col("status").eq(lit("collapse_other"))),
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
            .otherwise(col("status"))
            .alias("status"),
        )
        .with_columns([
            // misjoin
            // Regions with either:
            // * Zero coverage. Might be scaffold or misjoined contig. Without reference, not known.
            // * A dip in coverage with a high het ratio.
            when(col("first").eq(lit(0)).or(col("first_peak").eq(lit("low"))))
                .then(lit("misjoin"))
                // false_dupe
                // Region with half of the expected coverage, n-zscores less than the global mean, and zero-mapq due to multi-mapping.
                // Either a duplicated contig or duplicated region.
                .when(
                    col("first")
                        .lt_eq(lit(median_cov / 2))
                        // // Needs to scale with coverage.
                        // .and(col("first_all_zscore"))
                        .and(col("mapq").eq(lit(0))),
                )
                .then(lit("false_dupe"))
                .otherwise(col("status"))
                .alias("status"),
        ])
        .with_column(lit(ctg).alias("chrom"))
        .select([
            col("chrom"),
            col("pos"),
            col("first"),
            col("second"),
            col("first_zscore"),
            col("second_zscore"),
            col("mapq"),
            col("status"),
        ])
        .collect()?;

    // Construct intervals.
    // Store [st,end,type,cov]
    let df_itvs = df_pileup
        .select(["pos", "first", "second", "status"])?
        .lazy()
        .with_column(col("status").rle_id().alias("group"))
        .group_by([col("group")])
        .agg([
            col("pos").min().alias("st"),
            col("pos").max().alias("end") + lit(1),
            (col("first") + col("second"))
                .mean()
                .alias("cov")
                .cast(DataType::UInt8),
            col("status").first(),
        ])
        .drop([col("group")])
        .collect()?;
    
    Ok(
        (df_itvs, cfg.general.store_coverage.then_some(df_pileup))
    )
}

/// Detect misasemblies from alignment read coverage using per-base read coverage.
///
/// # Arguments
/// * `bamfile`: Input BAM file path. Should be indexed and filtered to only primary alignments (?).  
/// * `itv`: Interval to check.
/// * `cfg`: Peak-calling configuration.
/// * `avg_cov`: Average coverage of assembly. Used only for classifying false-duplications. Defaults to contig coverage.
///
/// # Returns
/// * [`NucFlagResult`]
pub fn nucflag(
    bamfile: impl AsRef<Path>,
    itv: &Interval<String>,
    cfg: Config,
    avg_cov: Option<u64>,
) -> eyre::Result<NucFlagResult> {
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.first.try_into()?, itv.last.try_into()?);
    let mut bam = bam::io::indexed_reader::Builder::default().build_from_path(&bamfile)?;
    let pileup = pileup(&mut bam, itv)?;

    // Scale thresholds by contig depth. Regions with fewer mapped reads get stricter parameters.
    let df_raw_pileup = merge_pileup_info(pileup.pileups, st, end)?;
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");

    let median_first: u64 = df_raw_pileup
        .column("first")?
        .median_reduce()?
        .value()
        .try_extract()?;
    let median_second: u64 = df_raw_pileup
        .column("second")?
        .median_reduce()?
        .value()
        .try_extract()?;
    let median_cov = avg_cov.unwrap_or(median_first + median_second);
    log::debug!("Average coverage for {ctg}:{st}-{end}: {median_cov}");

    let lf_first_peaks = find_peaks(
        df_raw_pileup.select(["pos", "first"])?,
        cfg.first.n_zscores_low,
        cfg.first.n_zscores_high,
        None,
    )?;
    let lf_second_peaks = find_peaks(
        df_raw_pileup.select(["pos", "second"])?,
        cfg.second.n_zscores_high,
        cfg.second.n_zscores_high,
        Some(cfg.second.min_perc),
    )?;

    let lf_pileup = lf_first_peaks
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
        .with_column(
            // Calculate het ratio.
            (col("second").cast(DataType::Float32)
                / (col("first").cast(DataType::Float32) + col("second").cast(DataType::Float32)))
            .alias("het_ratio"),
        );

    std::mem::drop(df_raw_pileup);

    let (df_itvs, df_pileup) = classify_peaks(lf_pileup, &ctg, &cfg, median_cov)?;

    let bp_filter: i32 = cfg.general.min_bp.try_into()?;
    let bp_merge: i32 = cfg.general.bp_merge.try_into()?;

    // Then merge and filter.
    log::info!("Merging intervals in {ctg}:{st}-{end}.");
    let df_itvs_final =
        merge_misassemblies(df_itvs, bp_merge, bp_filter, cfg.general.merge_across_type)?
            .with_columns([
                lit(ctg.clone()).alias("chrom"),
                col("st").alias("thickStart"),
                col("end").alias("thickEnd"),
                lit("+").alias("strand"),
                // Convert statuses into colors.
                col("status").map(
                    |statuses| {
                        Ok(Some(Column::new(
                            "itemRgb".into(),
                            statuses
                                .str()?
                                .iter()
                                .flatten()
                                .map(|s| MisassemblyType::from_str(s).unwrap().item_rgb())
                                .collect::<Vec<&str>>(),
                        )))
                    },
                    SpecialEq::same_type(),
                ).alias("itemRgb"),
            ])
            .rename(
                ["st", "end", "status", "cov"],
                ["chromStart", "chromEnd", "name", "score"],
                true,
            )
            .select([
                col("chrom"),
                col("chromStart"),
                col("chromEnd"),
                col("name"),
                col("score"),
                col("strand"),
                col("thickStart"),
                col("thickEnd"),
                col("itemRgb"),
            ])
            .collect()?;

    let (n_misassemblies, _) = df_itvs_final
        .select(["name"])?
        .lazy()
        .filter(col("name").neq(lit("good")))
        .collect()?
        .shape();

    log::info!("Detected {n_misassemblies} misassemblies for {ctg}:{st}-{end}.",);
    Ok(NucFlagResult {
        cov: df_pileup,
        regions: df_itvs_final,
    })
}
