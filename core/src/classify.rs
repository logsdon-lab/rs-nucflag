use core::str;
use std::{collections::HashMap, path::Path, str::FromStr};

use crate::{
    binning::group_pileup_by_ani,
    config::{Config, MinimumSizeConfig},
    intervals::merge_intervals,
    misassembly::MisassemblyType,
    peak::find_peaks,
    pileup::{AlignmentFile, PileupInfo},
};
use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use itertools::{multizip, Itertools};
use polars::prelude::*;
use rayon::prelude::*;

// TODO: Write macro to simplify.
fn merge_pileup_info(
    pileup: Vec<PileupInfo>,
    st: u64,
    end: u64,
    cfg: &Config,
) -> eyre::Result<DataFrame> {
    let (
        mut cov_cnts,
        mut mismatch_cnts,
        mut mapq_cnts,
        mut supp_cnts,
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
    for p in pileup.into_iter() {
        cov_cnts.push(p.n_cov);
        mismatch_cnts.push(p.n_mismatch);
        mapq_cnts.push(p.median_mapq().unwrap_or(0));
        supp_cnts.push(p.n_supp);
        indel_cnts.push(p.n_indel);
        softclip_cnts.push(p.n_softclip);
    }
    
    let mut lf = DataFrame::new(vec![
        Column::new("pos".into(), st..end + 1),
        Column::new("cov".into(), cov_cnts),
        Column::new("mismatch".into(), mismatch_cnts),
        Column::new("mapq".into(), mapq_cnts),
        Column::new("supp".into(), supp_cnts),
        Column::new("indel".into(), indel_cnts),
        Column::new("softclip".into(), softclip_cnts),
    ])?.lazy();

    for (colname, window_size) in [
        ("cov", cfg.cov.rolling_mean_window),
        ("mismatch", cfg.mismatch.rolling_mean_window),
        ("indel", cfg.indel.rolling_mean_window),
    ] {
        if let Some(window_size) = window_size {
            lf = lf
            .with_column(col(colname).rolling_mean(RollingOptionsFixedWindow {
                window_size,
                center: true,
                ..Default::default()
            }))
        };
    }
    Ok(lf.collect()?)
}

fn merge_misassemblies(
    df_itvs: DataFrame,
    bp_merge: i32,
    cfg_min_size: &MinimumSizeConfig,
    merge_across_type: bool,
) -> eyre::Result<LazyFrame> {
    let itvs_all: Vec<(u64, u64, u8, String)> = multizip((
        df_itvs.column("st")?.u64()?.iter().flatten(),
        df_itvs.column("end")?.u64()?.iter().flatten(),
        df_itvs.column("cov")?.u8()?.iter().flatten(),
        df_itvs
            .column("status")?
            .str()?
            .iter()
            .flatten()
            .map(String::from),
    ))
    .collect();

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
        acc.extend(merge_intervals(
            grps.into_iter()
                .map(|(st, end, cov, status)| Interval::new(st as i32, end as i32, (status, cov))),
            bp_merge,
            |_, _| true,
            // Average coverage.
            |itv_1, itv_2| (itv_1.metadata.0, (itv_1.metadata.1 + itv_2.metadata.1) / 2),
            |itv: Interval<(&str, u8)>| Interval::new(itv.first, itv.last, (grp, itv.metadata.1)),
        ));
        acc
    });

    // Merge overlapping intervals OVER status type choosing largest misassembly type.
    let final_misasm_itvs: COITree<(&str, u8), usize> = COITree::new(&if merge_across_type {
        merge_intervals(
            merged_misasm_itvs.into_iter(),
            bp_merge,
            |_, _| true,
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
    });

    let thr_minimum_sizes: HashMap<&str, u64> = HashMap::from_iter([
        ("good", u64::MAX),
        ("collapse", cfg_min_size.collapse.try_into()?),
        ("false_dupe", cfg_min_size.false_dupe.try_into()?),
        ("indel", cfg_min_size.indel.try_into()?),
        ("low_quality", cfg_min_size.low_quality.try_into()?),
        ("misjoin", cfg_min_size.misjoin.try_into()?),
        ("softclip", cfg_min_size.softclip.try_into()?),
    ]);

    // Convert good columns to misassembly types.
    let mut sts = Vec::with_capacity(itvs_all.len());
    let mut ends = Vec::with_capacity(itvs_all.len());
    let mut covs = Vec::with_capacity(itvs_all.len());
    let mut statuses = Vec::with_capacity(itvs_all.len());
    let mut minimum_sizes = Vec::with_capacity(itvs_all.len());
    for (st, end, cov, status) in itvs_all {
        let st = st.try_into()?;
        let end = end.try_into()?;
        let mut largest_ovl: Option<(&str, i32)> = None;
        final_misasm_itvs.query(st, end, |ovl_itv| {
            // Fully contained.
            let contains_itv = ovl_itv.first <= st && ovl_itv.last >= end;
            match (largest_ovl, contains_itv) {
                (None, true) => largest_ovl = Some((ovl_itv.metadata.0, ovl_itv.len())),
                (Some((_, other_itv_len)), true) => {
                    // Take larger overlap as status
                    if ovl_itv.len() > other_itv_len {
                        largest_ovl = Some((ovl_itv.metadata.0, ovl_itv.len()))
                    }
                }
                _ => (),
            }
        });

        sts.push(st);
        ends.push(end);
        covs.push(cov);

        let Some((new_status, _)) = largest_ovl else {
            minimum_sizes.push(thr_minimum_sizes[status.as_str()]);
            statuses.push(status);
            continue;
        };

        minimum_sizes.push(thr_minimum_sizes[new_status]);
        statuses.push(new_status.to_owned());
    }
    let df_itvs_all = DataFrame::new(vec![
        Column::new("st".into(), sts),
        Column::new("end".into(), ends),
        Column::new("cov".into(), covs),
        Column::new("status".into(), statuses),
        Column::new("thr_size".into(), minimum_sizes),
    ])?;

    Ok(df_itvs_all
        .lazy()
        // First reduce interval groups to min/max.
        .with_column(col("status").rle_id().alias("group"))
        .group_by(["group"])
        .agg([
            col("st").min(),
            col("end").max(),
            col("cov").median(),
            col("status").first(),
            col("thr_size").first(),
        ])
        // Reset intervals smaller than filter
        .with_column(
            when((col("end") - col("st")).lt_eq(col("thr_size")))
                .then(lit("good"))
                .otherwise(col("status"))
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
    median_cov: u64,
) -> eyre::Result<(DataFrame, Option<DataFrame>)> {
    let lf_pileup = lf_pileup
        .with_column(
            // indel
            // Region with insertion or deletion or soft clip that has high indel ratio and has a peak.
            when(
                (col("indel").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                    .gt_eq(lit(cfg.indel.ratio_indel))
                    .and(col("indel_peak").eq(lit("high"))),
            )
            .then(lit("indel"))
            .when(
                (col("softclip").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                    .gt_eq(lit(cfg.softclip.ratio_softclip))
                    .and(col("softclip_peak").eq(lit("high"))),
            )
            .then(lit("softclip"))
            .otherwise(lit("good"))
            .alias("status"),
        )
        .with_column(
            // collapse
            // Regions with at double the coverage and high indel peak.
            when(
                col("cov_peak")
                    .eq(lit("high"))
                    .and(col("cov").gt_eq(lit(median_cov / 2)))
                    .and(col("indel_peak").eq(lit("high"))),
            )
            .then(lit("collapse"))
            // misjoin
            // Regions with zero coverage or a dip in coverage with a indels or softclipping accounting for majority of the dip.
            // Might be scaffold or misjoined contig. Without reference, not known.
            .when(
                col("cov").eq(lit(0)).or(col("cov_peak").eq(lit("low")).and(
                    col("status")
                        .eq(lit("indel"))
                        .or(col("status").eq(lit("softclip"))),
                )),
            )
            .then(lit("misjoin"))
            // false_dupe
            // Region with half of the expected coverage and zero-mapq due to multi-mapping.
            // Either a duplicated contig, duplicated region, or an SV (large insertion of repetive region).
            .when(
                col("cov")
                    .lt_eq(lit(median_cov / 2))
                    .and(col("mapq").eq(lit(0))),
            )
            .then(lit("false_dupe"))
            .otherwise(col("status"))
            .alias("status"),
        )
        .with_column(
            // low_quality
            // Regions with high mismatch peak and het ratio.
            when(
                col("mismatch_ratio")
                    .gt_eq(lit(cfg.mismatch.ratio_het))
                    .and(col("mismatch_peak").eq(lit("high"))),
            )
            .then(lit("low_quality"))
            .otherwise(col("status"))
            .alias("status"),
        );

    let df_pileup = lf_pileup
        .with_column(lit(ctg).alias("chrom"))
        .select([
            col("chrom"),
            col("pos"),
            col("cov"),
            col("mismatch"),
            col("mapq"),
            col("status"),
            col("indel"),
            col("supp"),
            col("softclip"),
            col("cov_zscore"),
            col("mismatch_zscore"),
            col("indel_zscore"),
            col("softclip_zscore"),
            col("bin"),
        ])
        .collect()?;

    // Construct intervals.
    // Store [st,end,type,cov]
    let df_itvs = df_pileup
        .select(["pos", "cov", "mismatch", "status"])?
        .lazy()
        .with_column(col("status").rle_id().alias("group"))
        .group_by([col("group")])
        .agg([
            col("pos").min().alias("st"),
            col("pos").max().alias("end") + lit(1),
            col("cov").mean().cast(DataType::UInt8),
            col("status").first(),
        ])
        .drop([col("group")])
        .collect()?;

    Ok((df_itvs, cfg.general.store_pileup.then_some(df_pileup)))
}

fn nucflag_grp(
    df_pileup: DataFrame,
    cfg: &Config,
    ctg: &str,
) -> eyre::Result<(DataFrame, Option<DataFrame>)> {
    // Calculate est coverage for region or use provided.
    let est_median_cov: u64 = df_pileup
        .column("cov")?
        .median_reduce()?
        .value()
        .try_extract()?;
    let median_cov = cfg.cov.baseline.unwrap_or(est_median_cov);

    //  Detect dips and peaks in coverage.
    let lf_cov_peaks = find_peaks(
        df_pileup.select(["pos", "cov"])?,
        cfg.cov.n_zscores_low,
        cfg.cov.n_zscores_high,
    )?;

    // Detect indel peaks.
    let lf_indel_peaks = find_peaks(
        df_pileup.select(["pos", "indel"])?,
        // Don't care about dips in indels.
        cfg.indel.n_zscores_high,
        cfg.indel.n_zscores_high,
    )?;
    let lf_cov_peaks = lf_cov_peaks.join(
        lf_indel_peaks,
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    // Call peaks in mismatch-base signal.
    let lf_mismatch_peaks = find_peaks(
        df_pileup.select(["pos", "mismatch"])?,
        cfg.mismatch.n_zscores_high,
        cfg.mismatch.n_zscores_high,
    )?;
    let lf_cov_peaks =lf_cov_peaks.join(
        lf_mismatch_peaks,
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    // Detect indel peaks.
    let lf_softclip_peaks = find_peaks(
        df_pileup.select(["pos", "softclip"])?,
        // Don't care about dips in indels.
        cfg.softclip.n_zscores_high,
        cfg.softclip.n_zscores_high,
    )?;
    let lf_cov_peaks = lf_cov_peaks.join(
        lf_softclip_peaks,
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    // Add mapq
    let lf_pileup = lf_cov_peaks.join(
        df_pileup.select(["pos", "mapq"])?.lazy(),
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    ).with_column(
        // Calculate het ratio.
        (col("mismatch").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
            .alias("mismatch_ratio"),
    );


    // Add supplementary and indel counts.
    let lf_pileup = lf_pileup.join(
        df_pileup
            .select(["pos", "supp", "indel", "softclip", "bin"])?
            .lazy(),
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    std::mem::drop(df_pileup);

    classify_peaks(lf_pileup, ctg, cfg, median_cov)
}

/// Detect misasemblies from alignment read coverage using per-base read coverage.
///
/// # Arguments
/// * `aln`: Input BAM/CRAM file path. Should be indexed.  
/// * `fasta`: Input BAM file path. Required with CRAM. Also used for region binning.  
/// * `itv`: Interval to check.
/// * `cfg`: Peak-calling configuration. See [`Preset`] for configuration based on sequencing data type.
///
/// # Returns
/// * [`NucFlagResult`]
pub fn nucflag(
    aln: impl AsRef<Path>,
    fasta: Option<impl AsRef<Path> + Clone>,
    itv: &Interval<String>,
    cfg: Config,
) -> eyre::Result<NucFlagResult> {
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.first.try_into()?, itv.last.try_into()?);

    let mut aln = AlignmentFile::new(aln, fasta.clone())?;
    let pileup = aln.pileup(itv)?;

    let df_raw_pileup = merge_pileup_info(pileup.pileups, st, end, &cfg)?;
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");

    let df_pileup_groups = if let (Some(fasta), Some(cfg_grp_by_ani)) = (fasta, &cfg.group_by_ani) {
        group_pileup_by_ani(df_raw_pileup, fasta, itv, cfg_grp_by_ani)?
            .partition_by(["bin"], true)?
    } else {
        vec![df_raw_pileup
            .lazy()
            .with_column(lit(0).alias("bin"))
            .collect()?]
    };

    let (dfs_itvs, dfs_pileup): (Vec<LazyFrame>, Vec<Option<LazyFrame>>) = df_pileup_groups
        .into_par_iter()
        .map(|df_pileup_grp| {
            let (df_itv, df_pileup) = nucflag_grp(df_pileup_grp, &cfg, &ctg).unwrap();
            (df_itv.lazy(), df_pileup.map(|df| df.lazy()))
        })
        .unzip();
    let df_itvs = concat(dfs_itvs, Default::default())?
        .sort(["st"], Default::default())
        .collect()?;

    let df_pileup = if cfg.general.store_pileup {
        Some(
            concat(
                dfs_pileup.into_iter().flatten().collect::<Vec<LazyFrame>>(),
                Default::default(),
            )?
            .sort(["chrom", "pos"], Default::default())
            .collect()?,
        )
    } else {
        None
    };

    let bp_merge: i32 = cfg.general.bp_merge.try_into()?;

    // Then merge and filter.
    log::info!("Merging intervals in {ctg}:{st}-{end}.");
    let df_itvs_final = merge_misassemblies(
        df_itvs,
        bp_merge,
        &cfg.minimum_size.unwrap_or_default(),
        cfg.general.merge_across_type,
    )?
    .with_columns([
        lit(ctg.clone()).alias("chrom"),
        col("st").alias("thickStart"),
        col("end").alias("thickEnd"),
        lit("+").alias("strand"),
        // Convert statuses into colors.
        col("status")
            .map(
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
            )
            .alias("itemRgb"),
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
