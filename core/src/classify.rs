use core::str;
use std::{collections::HashMap, fmt::Debug, fs::File, io::BufReader, path::Path, str::FromStr};

use crate::{
    binning::BinStats,
    config::Config,
    intervals::{merge_intervals, overlap_length, subtract_intervals},
    misassembly::MisassemblyType,
    repeats::{detect_largest_repeat, Repeat},
};
use coitrees::{COITree, Interval, IntervalTree};
use eyre::Context;
use itertools::{multizip, Itertools};
use noodles::{
    core::{Position, Region},
    fasta,
};
use polars::{frame::row::Row, prelude::*};

fn split_at_ignored_intervals<'a>(
    st: i32,
    end: i32,
    status: &'a MisassemblyType,
    itree_ignore_itvs: &COITree<String, usize>,
) -> Option<Vec<Interval<&'a MisassemblyType>>> {
    // Trim interval by ignored intervals.
    let mut all_ovls = vec![];
    itree_ignore_itvs.query(st, end, |itv| {
        all_ovls.push(Interval::new(itv.first, itv.last, &MisassemblyType::Null));
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

    // Merge overlapping misassembly intervals OVER status type choosing largest misassembly type.
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
    let final_misasm_itvs: COITree<(MisassemblyType, u32), usize> =
        COITree::new(itvs_misasm.iter());
    let thr_minimum_sizes: HashMap<MisassemblyType, u64> = (&cfg_min_size).try_into()?;

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

    // Convert good intervals overlapping misassembly types.
    // Detect repeats.
    // Remove ignored intervals.
    let mut reclassified_itvs_all: Vec<Interval<(MisassemblyType, u32, u64)>> =
        Vec::with_capacity(itvs_all.len());
    for (st, end, cov, status, bin) in itvs_all {
        let st = st.try_into()?;
        let end = end.try_into()?;
        let len = (end - st) as f32;
        let mut largest_ovl: Option<MisassemblyType> = None;
        let mtype = MisassemblyType::from_str(status)?;
        final_misasm_itvs.query(st, end, |ovl_itv| {
            let ovl_prop = overlap_length(st, end, ovl_itv.first, ovl_itv.last) as f32 / len;
            if ovl_prop < 0.5 {
                return;
            }
            match largest_ovl {
                None => largest_ovl = Some(ovl_itv.metadata.0),
                Some(other_itv_mtype) => {
                    // Take larger overlap as status
                    if other_itv_mtype.gt(&mtype) {
                        largest_ovl = Some(other_itv_mtype)
                    }
                }
            }
        });
        let status = largest_ovl
            .filter(|ovl_mtype| ovl_mtype.gt(&mtype))
            .unwrap_or_else(|| MisassemblyType::from_str(status).unwrap());

        // Detect scaffold/homopolymer/simple repeat and replace type.
        let status = if let (Some(reader), Some(cfg_rpt)) = (
            fasta_reader.as_mut(),
            // Must have repeat config and the current status must be in types to check.
            cfg.repeat
                .as_ref()
                .map(|cfg_rpt| (mtype, cfg_rpt))
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
                .and_then(|rpt| {
                    log::debug!("Detected repeat at {ctg}:{st}-{end}: {rpt:?}");
                    // If any number of N's is scaffold.
                    if rpt.repeat == Repeat::Scaffold {
                        Some(MisassemblyType::Repeat(rpt.repeat))
                    } else {
                        (rpt.prop > cfg_rpt.ratio_repeat)
                            .then_some(MisassemblyType::Repeat(rpt.repeat))
                    }
                })
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
                reclassified_itvs_all.push(Interval::new(itv.first, itv.last, (status, cov, bin)));
            }
            continue;
        }

        // Otherwise, add misassembly.
        reclassified_itvs_all.push(Interval::new(st, end, (status, cov, bin)));
    }
    // Get minimum and maximum positions of sorted, grouped intervals.
    // Filter collapses based on calculate derivative within misassembled regions
    let mut minmax_reclassified_itvs_all = vec![];
    for ((is_mergeable, bin), group_elements) in &reclassified_itvs_all
        .into_iter()
        .sorted_by(|a, b| a.first.cmp(&b.first))
        .chunk_by(|a| (a.metadata.0.is_mergeable(), a.metadata.2))
    {
        if is_mergeable {
            let (mut agg_st, mut agg_end, mut mean_cov) = (i32::MAX, 0, 0);
            let mut agg_status = MisassemblyType::Null;
            let mut num_elems = 0;
            let elems: Vec<(i32, i32, MisassemblyType, u32)> = group_elements
                .map(|itv| (itv.first, itv.last, itv.metadata.0, itv.metadata.1))
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect();
            // Get min max of region.
            for (st, end, status, cov) in elems.iter() {
                agg_st = std::cmp::min(*st, agg_st);
                agg_end = std::cmp::max(*end, agg_end);
                agg_status = std::cmp::max(agg_status, *status);
                mean_cov += *cov;
                num_elems += 1;
            }
            mean_cov /= num_elems;

            // Change in coverage is greater than 1 stdev. Indicates that transition and should be ignored.
            // TODO: Might also apply to false dupe and misjoin. Check.
            let bin_stats = &bin_stats[&bin];
            if agg_status == MisassemblyType::Collapse {
                // Get runs of same signed changes.
                let min_max_change = elems
                    .windows(2)
                    .map(|a| a[1].3 as f32 - a[0].3 as f32)
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
                    agg_status = MisassemblyType::Null;
                }
            }
            minmax_reclassified_itvs_all.push(Interval::new(
                agg_st,
                agg_end,
                (agg_status, mean_cov),
            ));
        } else {
            minmax_reclassified_itvs_all.extend(
                group_elements.into_iter().map(|itv| {
                    Interval::new(itv.first, itv.last, (itv.metadata.0, itv.metadata.1))
                }),
            );
        }
    }

    // Remove intervals not within minimum sizes after merging.
    let minmax_reclassified_itvs_all = merge_intervals(
        minmax_reclassified_itvs_all.into_iter(),
        1,
        |a, b| a.metadata.0 == b.metadata.0,
        |a, b| (a.metadata.0, (a.metadata.1 + b.metadata.1) / 2),
        |a| {
            let mut status = a.metadata.0;
            // Remove misassemblies less than threshold size.
            let min_size = thr_minimum_sizes[&status];
            let length = (a.last - a.first) as u64;
            if length < min_size {
                status = MisassemblyType::Null;
            }
            Interval::new(a.first, a.last, (status, a.metadata.1))
        },
    );

    let minmax_reclassified_itvs_all: Vec<Row> = minmax_reclassified_itvs_all
        .into_iter()
        .map(|itv| {
            Row::new(vec![
                AnyValue::Int32(itv.first),
                AnyValue::Int32(itv.last),
                AnyValue::String(if itv.metadata.0 == MisassemblyType::Null {
                    "good"
                } else {
                    itv.metadata.0.into()
                }),
                AnyValue::UInt32(itv.metadata.1),
            ])
        })
        .collect();

    let df_itvs_all = DataFrame::from_rows_and_schema(
        &minmax_reclassified_itvs_all,
        &Schema::from_iter([
            ("st".into(), DataType::Int32),
            ("end".into(), DataType::Int32),
            ("status".into(), DataType::String),
            ("cov".into(), DataType::UInt32),
        ]),
    )?;

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
                .try_extract()
                .unwrap_or_default(),
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
