use core::str;
use std::collections::HashMap;

use crate::{
    config::{Config, MinimumSizeConfig},
    intervals::{merge_intervals, trim_coords},
    pileup::PileupInfo,
};
use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use itertools::{multizip, Itertools};
use polars::prelude::*;

pub(crate) fn merge_pileup_info(
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

pub(crate) fn merge_misassemblies(
    df_itvs: DataFrame,
    ignore_itvs: Option<&COITree<String, usize>>,
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
        let mut st = st.try_into()?;
        let mut end = end.try_into()?;
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

        let mut status = if let Some((new_status, _)) = largest_ovl {
            new_status.to_owned()
        } else {
            status
        };

        if let Some(itree_ignore_itvs) = ignore_itvs {
            let mut split_itv_1: Option<(i32, i32)> = None;
            let mut split_itv_2: Option<(i32, i32)> = None;
            let original_status = status.clone();
            // Trim interval by ignored intervals.
            itree_ignore_itvs.query(st, end, |itv| {
                trim_coords(
                    &mut st,
                    &mut end,
                    itv,
                    &mut status,
                    &mut split_itv_1,
                    &mut split_itv_2,
                )
            });

            // Add split intervals if any.
            if let (Some((st_itv_1, end_itv_1)), Some((st_itv_2, end_itv_2))) =
                (split_itv_1, split_itv_2)
            {
                sts.push(st_itv_1);
                ends.push(end_itv_1);
                covs.push(cov);
                sts.push(st_itv_2);
                ends.push(end_itv_2);
                covs.push(cov);

                minimum_sizes.push(thr_minimum_sizes[original_status.as_str()]);
                minimum_sizes.push(thr_minimum_sizes[original_status.as_str()]);
                statuses.push(original_status.clone());
                statuses.push(original_status);
            }
        }

        sts.push(st);
        ends.push(end);
        covs.push(cov);

        minimum_sizes.push(thr_minimum_sizes[status.as_str()]);
        statuses.push(status);
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

pub(crate) fn classify_peaks(
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
