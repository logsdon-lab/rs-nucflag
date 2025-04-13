use core::str;
use std::{path::Path, str::FromStr};

use crate::{
    config::Config, intervals::merge_overlapping_intervals, misassembly::MisassemblyType,
    peak::find_peaks, pileup::AlignmentFile, pileup::PileupInfo,
};
use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use itertools::{multizip, Itertools};
use noodles::{
    core::{Position, Region},
    fasta,
};
use polars::prelude::*;

// TODO: Write macro to simplify.
fn merge_pileup_info(
    pileup: Vec<PileupInfo>,
    st: u64,
    end: u64,
    cfg: &Config,
) -> eyre::Result<DataFrame> {
    let cols = if cfg.mismatch.is_some() {
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
        vec![
            Column::new("pos".into(), st..end + 1),
            Column::new("cov".into(), cov_cnts),
            Column::new("mismatch".into(), mismatch_cnts),
            Column::new("mapq".into(), mapq_cnts),
            Column::new("supp".into(), supp_cnts),
            Column::new("indel".into(), indel_cnts),
            Column::new("softclip".into(), softclip_cnts),
        ]
    } else {
        let (mut cov_cnts, mut mapq_cnts, mut supp_cnts, mut indel_cnts, mut softclip_cnts) = (
            Vec::with_capacity(pileup.len()),
            Vec::with_capacity(pileup.len()),
            Vec::with_capacity(pileup.len()),
            Vec::with_capacity(pileup.len()),
            Vec::with_capacity(pileup.len()),
        );
        for p in pileup.into_iter() {
            cov_cnts.push(p.n_cov);
            mapq_cnts.push(p.median_mapq().unwrap_or(0));
            supp_cnts.push(p.n_supp);
            indel_cnts.push(p.n_indel);
            softclip_cnts.push(p.n_softclip);
        }
        vec![
            Column::new("pos".into(), st..end + 1),
            Column::new("cov".into(), cov_cnts),
            Column::new("mapq".into(), mapq_cnts),
            Column::new("supp".into(), supp_cnts),
            Column::new("indel".into(), indel_cnts),
            Column::new("softclip".into(), softclip_cnts),
        ]
    };
    let df = DataFrame::new(cols)?;
    if let Some(window_size) = cfg.indel.rolling_mean_window {
        Ok(df
            .lazy()
            .with_column(col("indel").rolling_mean(RollingOptionsFixedWindow {
                window_size,
                center: true,
                ..Default::default()
            }))
            .collect()?)
    } else {
        Ok(df)
    }
}

fn merge_misassemblies(
    df_itvs: DataFrame,
    bp_merge: i32,
    bp_filter: i32,
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
    let final_misasm_itvs: COITree<(&str, u8), usize> = COITree::new(&if merge_across_type {
        merge_overlapping_intervals(
            merged_misasm_itvs
                .into_iter()
                .map(|itv| Interval::new(itv.first - bp_merge, itv.last + bp_merge, itv.metadata)),
            |itv_1, itv_2| {
                let largest_itv =
                    std::cmp::max_by(itv_1, itv_2, |itv_1, itv_2| itv_1.len().cmp(&itv_2.len()));
                (
                    largest_itv.metadata.0,
                    (itv_1.metadata.1 + itv_2.metadata.1) / 2,
                )
            },
            |itv: Interval<(&str, u8)>| {
                Interval::new(itv.first + bp_merge, itv.last - bp_merge, itv.metadata)
            },
        )
    } else {
        merged_misasm_itvs
    });

    // Convert good columns to misassembly types.
    let mut sts = Vec::with_capacity(itvs_all.len());
    let mut ends = Vec::with_capacity(itvs_all.len());
    let mut covs = Vec::with_capacity(itvs_all.len());
    let mut statuses = Vec::with_capacity(itvs_all.len());
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
            statuses.push(status);
            continue;
        };

        statuses.push(new_status.to_owned());
    }
    let df_itvs_all = DataFrame::new(vec![
        Column::new("st".into(), sts),
        Column::new("end".into(), ends),
        Column::new("cov".into(), covs),
        Column::new("status".into(), statuses),
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
            // Regions with at high coverage and high indel peak.
            when(
                col("cov_peak")
                    .eq(lit("high"))
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
            // Region with half of the expected coverage, n-zscores less than the global mean, and zero-mapq due to multi-mapping.
            // Either a duplicated contig, duplicated region, or an SV (large insertion of repetive region).
            .when(
                col("cov")
                    .lt_eq(lit(median_cov / 2))
                    .and(col("cov_zscore").lt_eq(lit(-cfg.cov.n_zscores_false_dupe)))
                    .and(col("mapq").eq(lit(0))),
            )
            .then(lit("false_dupe"))
            .otherwise(col("status"))
            .alias("status"),
        );

    let lf_pileup = if let Some(cfg_mismatch) = &cfg.mismatch {
        lf_pileup.with_column(
            // low_quality
            // Regions with high mismatch peak and het ratio.
            when(
                col("mismatch_ratio")
                    .gt_eq(lit(cfg_mismatch.ratio_het))
                    .and(col("mismatch_peak").eq(lit("high"))),
            )
            .then(lit("low_quality"))
            .otherwise(col("status"))
            .alias("status"),
        )
    } else {
        lf_pileup
    };

    let pileup_cols = if cfg.mismatch.is_some() {
        vec![
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
        ]
    } else {
        vec![
            col("chrom"),
            col("pos"),
            col("cov"),
            col("mapq"),
            col("status"),
            col("indel"),
            col("supp"),
            col("softclip"),
            col("cov_zscore"),
            col("indel_zscore"),
            col("softclip_zscore"),
        ]
    };
    let df_pileup = lf_pileup
        .with_column(lit(ctg).alias("chrom"))
        .select(pileup_cols)
        .collect()?;

    // Construct intervals.
    // Store [st,end,type,cov]
    let itv_cols = if cfg.mismatch.is_some() {
        vec!["pos", "cov", "mismatch", "status"]
    } else {
        vec!["pos", "cov", "status"]
    };
    let df_itvs = df_pileup
        .select(itv_cols)?
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

    Ok((df_itvs, cfg.general.store_coverage.then_some(df_pileup)))
}

/// Detect misasemblies from alignment read coverage using per-base read coverage.
///
/// # Arguments
/// * `aln`: Input BAM/CRAM file path. Should be indexed.  
/// * `fasta`: Input BAM file path. Required with CRAM. Also used for region binning.  
/// * `itv`: Interval to check.
/// * `cfg`: Peak-calling configuration.
/// * `avg_cov`: Average coverage of assembly. Used only for classifying false-duplications. Defaults to contig coverage.
///
/// # Returns
/// * [`NucFlagResult`]
pub fn nucflag(
    aln: impl AsRef<Path>,
    fasta: Option<impl AsRef<Path> + Clone>,
    itv: &Interval<String>,
    cfg: Config,
    avg_cov: Option<u64>,
) -> eyre::Result<NucFlagResult> {
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.first.try_into()?, itv.last.try_into()?);

    let mut aln = AlignmentFile::new(aln, fasta.clone())?;
    let pileup = aln.pileup(itv)?;

    let df_raw_pileup = merge_pileup_info(pileup.pileups, st, end, &cfg)?;
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");

    let df_raw_pileup = if let Some(fasta) = fasta {
        let mut reader_fasta =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;
        let position = Position::new(itv.first.try_into()?).unwrap()
            ..=Position::new(itv.last.try_into()?).unwrap();
        let region = Region::new(itv.metadata.clone(), position);
        let seq = reader_fasta.query(&region)?;
        eprintln!("Calculating self-identity for {ctg}:{st}-{end} to bin region.");
        // TODO: Adjust coordinates.
        let bed_ident = rs_moddotplot::compute_seq_self_identity(
            str::from_utf8(seq.sequence().as_ref())?,
            &itv.metadata,
            None,
        );
        println!("{bed_ident:?}");
        df_raw_pileup
    } else {
        df_raw_pileup
    };

    // Calculate est coverage for region or use provided.
    let est_median_cov: u64 = df_raw_pileup
        .column("cov")?
        .median_reduce()?
        .value()
        .try_extract()?;
    let median_cov = avg_cov.unwrap_or(est_median_cov);
    log::debug!("Average coverage for {ctg}:{st}-{end}: {median_cov}");

    //  Detect dips and peaks in coverage.
    let lf_cov_peaks = find_peaks(
        df_raw_pileup.select(["pos", "cov"])?,
        cfg.cov.n_zscores_low,
        cfg.cov.n_zscores_high,
    )?;

    // Detect indel peaks.
    let lf_indel_peaks = find_peaks(
        df_raw_pileup.select(["pos", "indel"])?,
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

    // Optionally, call peaks in mismatch-base signal.
    let lf_cov_peaks = if let Some(cfg_mismatch) = &cfg.mismatch {
        let lf_mismatch_peaks = find_peaks(
            df_raw_pileup.select(["pos", "mismatch"])?,
            cfg_mismatch.n_zscores_high,
            cfg_mismatch.n_zscores_high,
        )?;
        lf_cov_peaks.join(
            lf_mismatch_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
    } else {
        lf_cov_peaks
    };

    // Detect indel peaks.
    let lf_softclip_peaks = find_peaks(
        df_raw_pileup.select(["pos", "softclip"])?,
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
        df_raw_pileup.select(["pos", "mapq"])?.lazy(),
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    // Calculate het ratio.
    let lf_pileup = if cfg.mismatch.is_some() {
        lf_pileup.with_column(
            (col("mismatch").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                .alias("mismatch_ratio"),
        )
    } else {
        lf_pileup
    };

    // Add supplementary and indel counts.
    let lf_pileup = lf_pileup.join(
        df_raw_pileup
            .select(["pos", "supp", "indel", "softclip"])?
            .lazy(),
        [col("pos")],
        [col("pos")],
        JoinArgs::new(JoinType::Left),
    );

    std::mem::drop(df_raw_pileup);

    let (df_itvs, df_pileup) = classify_peaks(lf_pileup, &ctg, &cfg, median_cov)?;

    let bp_filter: i32 = cfg.general.bp_min.try_into()?;
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
