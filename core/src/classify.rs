use core::str;
use std::{collections::HashMap, fmt::Debug, fs::File, io::BufReader, path::Path, str::FromStr};

use crate::{
    binning::BinStats,
    config::Config,
    intervals::{merge_intervals, subtract_intervals},
    misassembly::MisassemblyType,
    repeats::detect_largest_repeat,
};
use coitrees::{COITree, Interval, IntervalTree};
use eyre::Context;
use itertools::{multizip, Itertools};
use noodles::{
    core::{Position, Region},
    fasta,
};
use polars::prelude::*;

fn split_at_ignored_intervals<'a>(
    st: i32,
    end: i32,
    status: &'a str,
    itree_ignore_itvs: &COITree<String, usize>,
) -> Option<Vec<Interval<&'a str>>> {
    // Trim interval by ignored intervals.
    let mut all_ovls = vec![];
    itree_ignore_itvs.query(st, end, |itv| {
        all_ovls.push(Interval::new(itv.first, itv.last, ""));
    });

    if all_ovls.is_empty() {
        return None;
    }

    let curr_itv = Interval::new(st, end, status);
    let new_itvs = subtract_intervals(curr_itv, all_ovls.into_iter());

    // If not equal to initial interval, nothing overlaps. Allow through.
    if new_itvs
        .first()
        .is_some_and(|i| i.first == curr_itv.first && i.last == curr_itv.last)
    {
        None
    } else {
        Some(new_itvs)
    }
}

pub(crate) fn merge_misassemblies(
    df_itvs: DataFrame,
    bin_stats: HashMap<u64, BinStats>,
    ctg: &str,
    fasta: Option<impl AsRef<Path> + Debug>,
    ignore_itvs: Option<&COITree<String, usize>>,
    cfg: Config,
) -> eyre::Result<LazyFrame> {
    let bp_merge = cfg.general.bp_merge.try_into()?;
    let cfg_min_size = cfg.minimum_size.unwrap_or_default();
    let itvs_all: Vec<(u64, u64, u32, &str, u64)> = multizip((
        df_itvs.column("st")?.u64()?.iter().flatten(),
        df_itvs.column("end")?.u64()?.iter().flatten(),
        df_itvs.column("cov")?.u32()?.iter().flatten(),
        df_itvs.column("status")?.str()?.iter().flatten(),
        df_itvs.column("bin")?.u64()?.iter().flatten(),
    ))
    .collect();

    let df_misasm_itvs = df_itvs
        .clone()
        .lazy()
        .filter(col("status").neq(lit("good")))
        .collect()?;

    let itvs_misasm = merge_intervals(
        multizip((
            df_misasm_itvs.column("st")?.u64()?.iter().flatten(),
            df_misasm_itvs.column("end")?.u64()?.iter().flatten(),
            df_misasm_itvs.column("cov")?.u32()?.iter().flatten(),
            df_misasm_itvs.column("status")?.str()?.iter().flatten(),
        ))
        .map(|(st, end, cov, status)| {
            Interval::new(
                st as i32,
                end as i32,
                (MisassemblyType::from_str(status).unwrap(), cov),
            )
        }),
        bp_merge,
        |_, _| true,
        |itv_1, itv_2| {
            let largest_itv = std::cmp::max_by(itv_1, itv_2, |itv_1, itv_2| {
                itv_1.metadata.0.cmp(&itv_2.metadata.0)
            });
            (
                largest_itv.metadata.0,
                (itv_1.metadata.1 + itv_2.metadata.1) / 2,
            )
        },
        |itv| itv,
    );
    // Merge overlapping intervals OVER status type choosing largest misassembly type.
    let final_misasm_itvs: COITree<(MisassemblyType, u32), usize> =
        COITree::new(itvs_misasm.iter());
    let thr_minimum_sizes: HashMap<&str, u64> = (&cfg_min_size).try_into()?;

    let mut fasta_reader = if let Some(fasta) = fasta {
        log::info!("Reading indexed {fasta:?} for {ctg} to detect further classify misassemblies by repeat.");
        let mut fai = fasta.as_ref().as_os_str().to_owned();
        fai.push(".fai");
        let fai =
            fasta::fai::read(fai).with_context(|| format!("Fasta file {fasta:?} not indexed."))?;
        let fa = BufReader::new(File::open(fasta)?);
        Some(fasta::IndexedReader::new(fa, fai))
    } else {
        None
    };

    // Convert good columns to misassembly types.
    let mut sts = Vec::with_capacity(itvs_all.len());
    let mut ends = Vec::with_capacity(itvs_all.len());
    let mut covs = Vec::with_capacity(itvs_all.len());
    let mut statuses = Vec::with_capacity(itvs_all.len());
    let mut bins = Vec::with_capacity(itvs_all.len());
    for (st, end, cov, status, bin) in itvs_all {
        let st = st.try_into()?;
        let end = end.try_into()?;
        let mut largest_ovl: Option<(MisassemblyType, i32)> = None;
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

        let status = if let Some((new_status, _)) = largest_ovl {
            Into::<&'static str>::into(new_status).to_owned()
        } else {
            status.to_owned()
        };

        // Detect scaffold/homopolymer/simple repeat and replace type.
        let status = if let (Some(reader), Some(cfg_rpt)) = (
            fasta_reader.as_mut(),
            // Must have repeat config and the current status must be in types to check.
            MisassemblyType::from_str(&status)
                .map(|mtype| cfg.repeat.as_ref().map(|cfg_rpt| (mtype, cfg_rpt)))?
                .and_then(|(mtype, cfg_rpt)| {
                    cfg_rpt.check_types.contains(&mtype).then_some(cfg_rpt)
                }),
        ) {
            // Add extended region.
            let st = st
                .saturating_sub(cfg_rpt.bp_extend.try_into()?)
                .clamp(1, i32::MAX)
                .try_into()?;
            let end = end
                .saturating_add(cfg_rpt.bp_extend.try_into()?)
                .try_into()?;
            let region = Region::new(
                ctg,
                Position::new(st).unwrap()..=Position::new(end).unwrap(),
            );
            let record = reader.query(&region)?;
            let seq = str::from_utf8(record.sequence().as_ref())?;
            detect_largest_repeat(seq)
                .and_then(|rpt| (rpt.prop > cfg_rpt.ratio_repeat).then(|| rpt.repeat.to_string()))
                .unwrap_or(status)
        } else {
            status
        };

        // This might not be the best approach, but it's the easiest :)
        // Ignoring during the pileup is better as it avoids even considering the region in calculations.
        // However, it complicates smoothing among other things.
        //
        // Split at ignored intervals if any overlap.
        if let Some(split_intervals) =
            ignore_itvs.and_then(|itree| split_at_ignored_intervals(st, end, &status, itree))
        {
            for itv in split_intervals {
                sts.push(itv.first);
                ends.push(itv.last);
                covs.push(cov);
                statuses.push(status.clone());
                bins.push(bin);
            }
            continue;
        }

        // Otherwise, add misassembly.
        sts.push(st);
        ends.push(end);
        covs.push(cov);
        statuses.push(status);
        bins.push(bin);
    }

    // Calculate derivative within misassembled regions to filter collapses
    let mut final_sts = Vec::with_capacity(sts.len());
    let mut final_ends = Vec::with_capacity(ends.len());
    let mut final_covs = Vec::with_capacity(ends.len());
    let mut final_statuses = Vec::with_capacity(statuses.len());
    let group_iter = multizip((sts, ends, covs, statuses, bins));
    for ((mut status, bin), group_elements) in &group_iter
        .sorted_by(|a, b| a.0.cmp(&b.0))
        .chunk_by(|a| (a.3.to_owned(), a.4))
    {
        let (mut agg_st, mut agg_end, mut mean_cov) = (i32::MAX, 0, 0);
        let mut num_elems = 0;
        let elems: Vec<(i32, i32, u32)> = group_elements
            .map(|(st, end, cov, _, _)| (st, end, cov))
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .collect();
        // Get min max of region.
        for (st, end, cov) in elems.iter() {
            agg_st = std::cmp::min(*st, agg_st);
            agg_end = std::cmp::max(*end, agg_end);
            mean_cov += *cov;
            num_elems += 1;
        }
        // Remove misassemblies less than threshold size.
        let min_size = thr_minimum_sizes[status.as_str()];
        let length = (agg_end - agg_st) as u64;
        if length < min_size {
            status.clear();
            status.push_str("good");
        }
        mean_cov /= num_elems;

        // Change in coverage is greater than 1 stdev. Indicates that transition and should be ignored.
        // TODO: Might also apply to false dupe and misjoin. Check.
        let bin_stats = &bin_stats[&bin];
        if status == "collapse" {
            // Get runs of same signed changes.
            let min_max_change = elems
                .windows(2)
                .map(|a| a[1].2 as f32 - a[0].2 as f32)
                .chunk_by(|a| a.is_sign_negative())
                .into_iter()
                .map(|(_, run)| run.sum::<f32>())
                .minmax();
            let stdev = bin_stats.stdev.abs();
            // eprintln!("Collapse: {ctg}:{agg_st}-{agg_end} with coverage changes ({min_max_change:?}) less than 1 stdev in bin {bin_stats:?}");
            if let Some((min, max)) = min_max_change
                .into_option()
                .filter(|(min, max)| (min.abs() < stdev.abs()) || (max.abs()) < stdev.abs())
            {
                log::debug!("Filtered out collapse: {ctg}:{agg_st}-{agg_end} with coverage changes ({min},{max}) less than 1 stdev in bin {bin_stats:?}");
                status.clear();
                status.push_str("good");
            }
        }
        final_sts.push(agg_st);
        final_ends.push(agg_end);
        final_covs.push(mean_cov);
        final_statuses.push(status);
    }

    let df_itvs_all = DataFrame::new(vec![
        Column::new("st".into(), final_sts),
        Column::new("end".into(), final_ends),
        Column::new("cov".into(), final_covs),
        Column::new("status".into(), final_statuses),
    ])?;

    // Reduce final interval groups to min/max.
    Ok(df_itvs_all
        .lazy()
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
    pub pileup: DataFrame,
}

pub(crate) fn classify_peaks(
    lf_pileup: LazyFrame,
    ctg: &str,
    cfg: &Config,
    median_cov: u32,
) -> eyre::Result<(DataFrame, DataFrame, BinStats)> {
    let thr_false_dupe = (cfg.cov.ratio_false_dupe * median_cov as f32).floor();
    let thr_collapse = (cfg.cov.ratio_collapse * median_cov as f32).floor();
    let thr_misjoin = (cfg.cov.ratio_misjoin * median_cov as f32).floor();

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
                    .and(col("cov").gt_eq(lit(thr_collapse)))
                    .and(col("indel_peak").eq(lit("high"))),
            )
            .then(lit("collapse"))
            // misjoin
            // Regions with zero coverage or a dip in coverage with a indels or softclipping accounting for majority of the dip.
            // Might be scaffold or misjoined contig. Without reference, not known.
            .when(
                col("cov")
                    .eq(lit(0))
                    .or(col("cov_peak")
                        .eq(lit("low"))
                        .and(col("cov").lt_eq(thr_misjoin)))
                    .or((col("cov") - col("indel")).lt_eq(thr_misjoin)),
            )
            .then(lit("misjoin"))
            // false_dupe
            // Region with half of the expected coverage and zero-mapq due to multi-mapping.
            // Either a duplicated contig, duplicated region, or an SV (large insertion of repetive region).
            .when(
                col("cov")
                    .lt_eq(lit(thr_false_dupe))
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
                (col("mismatch").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                    .gt_eq(lit(cfg.mismatch.ratio_het))
                    .and(col("mismatch_peak").eq(lit("high"))),
            )
            .then(lit("low_quality"))
            .otherwise(col("status"))
            .alias("status"),
        );

    let bin_stats = {
        // Apply a rolling median and stdev to get bin statistics.
        let rolling_opts = RollingOptionsFixedWindow {
            window_size: cfg.general.bp_merge,
            center: true,
            ..Default::default()
        };
        let df_bin_stats = lf_pileup
            .clone()
            // Only use good regions for bin stats.
            .filter(col("status").eq(lit("good")))
            .with_columns([
                col("cov")
                    .rolling_median(rolling_opts.clone())
                    .alias("cov_median"),
                col("cov").rolling_std(rolling_opts).alias("cov_stdev"),
            ])
            .select([col("bin"), col("cov_median"), col("cov_stdev")])
            .collect()?;
        BinStats {
            num: df_bin_stats
                .column("bin")?
                .u64()?
                .first()
                .unwrap_or_default(),
            median: df_bin_stats
                .column("cov_median")?
                .median_reduce()?
                .value()
                .try_extract()?,
            stdev: df_bin_stats
                .column("cov_stdev")?
                .median_reduce()?
                .value()
                .try_extract()?,
        }
    };
    /*
    // Removed cols to reduce memory consumption.
    col("cov_zscore"),
    col("mismatch_zscore"),
    col("indel_zscore"),
    col("softclip_zscore"),
    */
    let cols = [
        col("chrom"),
        col("pos"),
        col("cov"),
        col("status"),
        col("mismatch"),
        col("mapq"),
        col("indel"),
        col("softclip"),
        col("bin"),
    ];
    let df_pileup = lf_pileup
        .with_column(lit(ctg).alias("chrom"))
        .select(cols)
        .collect()?;

    // Construct intervals.
    // Store [st,end,type,cov]
    let df_itvs = df_pileup
        .select(["pos", "cov", "mismatch", "status", "bin"])?
        .lazy()
        .with_column(
            ((col("pos") - col("pos").shift_and_fill(1, 0))
                .lt_eq(1)
                .rle_id()
                + col("status").rle_id())
            .alias("group"),
        )
        .group_by([col("group")])
        .agg([
            col("pos").min().alias("st"),
            col("pos").max().alias("end") + lit(1),
            col("cov").mean().cast(DataType::UInt32),
            col("status").first(),
            col("bin").first(),
        ])
        .drop([col("group")])
        .collect()?;

    Ok((df_itvs, df_pileup, bin_stats))
}
