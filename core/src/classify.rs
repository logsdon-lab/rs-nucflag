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
use eyre::{bail, Context};
use itertools::{multizip, Itertools};
use noodles::{
    core::{Position, Region},
    fasta::{self, IndexedReader},
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

fn get_itree_above_median(
    lf_pileup: LazyFrame,
    median_cov: u32,
) -> eyre::Result<COITree<(), usize>> {
    let df_above_median_cov = lf_pileup
        .select([col("pos"), col("cov")])
        .with_columns([
            col("pos").min().alias("min_pos"),
            col("pos").max().alias("max_pos"),
            col("cov").gt_eq(lit(median_cov)).alias("above_median"),
        ])
        // +++----+++
        // 0001111222
        // 000----222
        // 012----789
        .with_column(col("above_median").rle_id())
        .filter(
            col("cov")
                .gt_eq(lit(median_cov))
                .over([col("above_median")]),
        )
        .group_by([col("above_median")])
        .agg([
            col("pos").min().alias("start"),
            col("pos").max().alias("end"),
            col("min_pos").first(),
            col("max_pos").first(),
        ])
        // Does the interval above the median contain the min or max position?
        .filter(
            col("min_pos")
                .gt_eq(col("start"))
                .and(col("min_pos").lt_eq(col("end")))
                .or(col("max_pos")
                    .gt_eq(col("start"))
                    .and(col("max_pos").lt_eq(col("end")))),
        )
        .select([col("start"), col("end")])
        .collect()?;
    let itvs_above_median_cov: Vec<Interval<()>> = df_above_median_cov
        .column("start")?
        .u64()?
        .iter()
        .flatten()
        .zip(df_above_median_cov.column("end")?.u64()?.iter().flatten())
        .map(|(st, end)| Interval::new(st as i32, end as i32, ()))
        .collect();
    Ok(COITree::new(&itvs_above_median_cov))
}

fn ignore_boundary_misassemblies(
    itvs: &mut [Interval<(MisassemblyType, u32, u64)>],
    ctg: &str,
    fasta: Option<IndexedReader<BufReader<File>>>,
    bin_stats: &HashMap<u64, BinStats>,
    default_boundary_positions: (i32, i32),
) {
    // Filter boundary misassemblies if below median coverage.
    // Handles cases in telomeres classified as misjoin/false_dupe/indel/repeats just because fewer reads.
    // * Also useful for specific regions like centromeres where we only care about the active array and don't mind misassemblies in pericentromere.
    let (ctg_st, ctg_end) = fasta
        .as_ref()
        .and_then(|fh| {
            let rec = fh.index().as_ref().iter().find(|rec| rec.name() == ctg)?;
            Some((0i32, (rec.length() + 1) as i32))
        })
        .unwrap_or(default_boundary_positions);

    // Keep going if below median in both directions.
    // With median of 8x
    // coverage  0 1 2 8  8  8 3 0
    // status  | x x x o ... o x x |
    // Each x is replace with a Null classification.
    let mut idx_st = 0;
    let mut idx_end = itvs.len() - 1;

    // Check for contig/queried region start or end.
    if itvs
        .first()
        .map(|itv| itv.first == ctg_st)
        .unwrap_or_default()
    {
        while let Some(itv) = itvs.get_mut(idx_st).filter(|itv| {
            let bin = &bin_stats[&itv.metadata.2];
            itv.metadata.1 < bin.median as u32
        }) {
            let og_mdata = itv.metadata;
            log::debug!("Filtered out {:?}: {ctg}:{}-{} at contig start with coverage below bin median {:?}", og_mdata.0, itv.first, itv.last, &bin_stats[&og_mdata.2]);
            *itv = Interval::new(
                itv.first,
                itv.last,
                (MisassemblyType::Null, og_mdata.1, og_mdata.2),
            );
            idx_st += 1
        }
    }

    if itvs
        .last()
        .map(|itv| itv.last == ctg_end)
        .unwrap_or_default()
    {
        while let Some(itv) = itvs.get_mut(idx_end).filter(|itv| {
            let bin = &bin_stats[&itv.metadata.2];
            itv.metadata.1 < bin.median as u32
        }) {
            let og_mdata = itv.metadata;
            log::debug!(
                "Filtered out {:?}: {ctg}:{}-{} on contig end with coverage below bin median {:?}",
                og_mdata.0,
                itv.first,
                itv.last,
                &bin_stats[&og_mdata.2]
            );
            *itv = Interval::new(
                itv.first,
                itv.last,
                (MisassemblyType::Null, og_mdata.1, og_mdata.2),
            );
            idx_end -= 1
        }
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

    let (Some(all_st), Some(all_end)) = (
        itvs_all.first().map(|itv| itv.0 as i32),
        itvs_all.last().map(|itv| itv.1 as i32),
    ) else {
        bail!("No intervals for {ctg}. Something is wrong.");
    };

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
        log::info!("Reading indexed {fasta:?} for {ctg} to classify misassemblies by repeat.");
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

        // Detect scaffold/homopolymer/repeat and replace type.
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
                        Some(MisassemblyType::RepeatError(rpt.repeat))
                    } else {
                        (rpt.prop > cfg_rpt.ratio_repeat)
                            .then_some(MisassemblyType::RepeatError(rpt.repeat))
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

    // Keep sorted.
    reclassified_itvs_all.sort_by(|a, b| a.first.cmp(&b.first));

    // Ignore boundary misassemblies.
    if cfg.general.ignore_boundaries {
        ignore_boundary_misassemblies(
            &mut reclassified_itvs_all,
            ctg,
            fasta_reader,
            &bin_stats,
            (all_st, all_end),
        );
    }

    // Get minimum and maximum positions of sorted, grouped intervals.
    // Filter collapses based on bin boundaries.
    let mut minmax_reclassified_itvs_all = vec![];
    for ((is_mergeable, bin), group_elements) in &reclassified_itvs_all
        .into_iter()
        .chunk_by(|a| (a.metadata.0.is_mergeable(), a.metadata.2))
    {
        if is_mergeable {
            let (mut agg_st, mut agg_end, mut mean_cov) = (i32::MAX, 0, 0);
            let mut agg_status = MisassemblyType::Null;
            let mut num_elems = 0;
            // Get min max of region.
            for (st, end, status, cov) in group_elements
                .map(|itv| (itv.first, itv.last, itv.metadata.0, itv.metadata.1))
                .sorted_by(|a, b| a.0.cmp(&b.0))
            {
                agg_st = std::cmp::min(st, agg_st);
                agg_end = std::cmp::max(end, agg_end);
                agg_status = std::cmp::max(agg_status, status);
                mean_cov += cov;
                num_elems += 1;
            }
            mean_cov /= num_elems;

            // At boundary of bin and is above median. Indicates that transition and should be ignored.
            //        v
            // 000000011111
            //    ____
            //  _/    \____
            // /           \
            let bin_stats = &bin_stats[&bin];
            if agg_status == MisassemblyType::Collapse
                && bin_stats
                    .itree_above_median
                    .query_count(agg_st, agg_end)
                    .ge(&1)
            {
                log::debug!("Filtered out {agg_status:?}: {ctg}:{agg_st}-{agg_end} above median coverage at bin transition ({bin_stats:?})");
                agg_status = MisassemblyType::Null;
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
            col("cov").median().cast(DataType::UInt32),
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
            // Regions with at double the coverage and dip in mapping quality or increase in mismatches/indels.
            when(
                col("cov_peak")
                    .eq(lit("high"))
                    .and(col("cov").gt_eq(lit(thr_collapse)))
                    .and(col("mapq_peak").eq(lit("low")))
                    .and(
                        col("mismatch_peak")
                            .eq(lit("high"))
                            .or(col("indel_peak").eq(lit("high"))),
                    ),
            )
            .then(lit("collapse"))
            // misjoin
            // Regions with zero coverage.
            .when(col("cov").eq(lit(0)))
            .then(lit("misjoin"))
            // false_dupe
            // Region with half of the expected coverage and a maximum mapq of zero due to multiple primary mappings.
            // Either a duplicated contig, duplicated region, or an SV (large insertion of repetive region).
            .when(
                col("cov")
                    .lt_eq(lit(thr_false_dupe))
                    .and(col("mapq_max").eq(lit(0))),
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
        let itree_above_median_cov = get_itree_above_median(lf_pileup.clone(), median_cov)?;
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
                .try_extract()
                .unwrap_or_default(),
            stdev: df_bin_stats
                .column("cov_stdev")?
                .median_reduce()?
                .value()
                .try_extract()
                .unwrap_or_default(),
            itree_above_median: itree_above_median_cov,
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
