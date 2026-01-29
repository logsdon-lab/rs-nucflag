use std::{collections::HashMap, fmt::Debug, path::Path, str::FromStr};

use crate::{
    binning::BinStats,
    config::Config,
    intervals::{merge_intervals, overlap_length, subtract_intervals},
    io::FastaHandle,
    misassembly::MisassemblyType,
    repeats::{detect_largest_repeat, Repeat},
};
use coitrees::{COITree, Interval, IntervalTree};
use eyre::bail;
use itertools::{multizip, Itertools};
use ordered_float::OrderedFloat;
use polars::lazy::dsl::max_horizontal;
use polars::{frame::row::Row, prelude::*};

#[derive(Debug, Clone)]
struct CallInfo {
    typ: MisassemblyType,
    cov: u32,
    bin: u32,
    zscore: OrderedFloat<f32>,
    af: OrderedFloat<f32>,
}

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
    itvs: &mut [Interval<CallInfo>],
    ctg: &str,
    fasta: Option<FastaHandle>,
    bin_stats: &HashMap<u32, BinStats>,
    default_boundary_positions: (i32, i32),
) {
    // Filter boundary misassemblies if below median coverage.
    // Handles cases in telomeres classified as misjoin/false_dup/indel/repeats just because fewer reads.
    // * Also useful for specific regions like centromeres where we only care about the active array and don't mind misassemblies in pericentromere.
    let (ctg_st, ctg_end) = fasta
        .as_ref()
        .and_then(|fh| {
            let rec = fh.fai.as_ref().iter().find(|rec| rec.name() == ctg)?;
            Some((0i32, (rec.length() + 1) as i32))
        })
        .unwrap_or(default_boundary_positions);

    // Keep going if below median - 1 stdev in both directions.
    // * Due to new merging rules, we cannot rely on presence of null to stop.
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
        // Keep removing while below median and not a good interval.
        while let Some(itv) = itvs.get_mut(idx_st).filter(|itv| {
            let bin = &bin_stats[&itv.metadata.bin];
            itv.metadata.cov < (bin.median - bin.stdev).clamp(0.0, f32::MAX) as u32
        }) {
            let og_mdata = &itv.metadata;
            log::debug!("Filtered out {:?}: {ctg}:{}-{} at contig start with coverage ({}) below bin median {:?}", og_mdata.typ, itv.first, itv.last, og_mdata.cov, &bin_stats[&og_mdata.bin]);
            *itv = Interval::new(
                itv.first,
                itv.last,
                CallInfo {
                    typ: MisassemblyType::Null,
                    cov: og_mdata.cov,
                    bin: og_mdata.bin,
                    zscore: og_mdata.zscore,
                    af: og_mdata.af,
                },
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
            let bin = &bin_stats[&itv.metadata.bin];
            itv.metadata.cov < (bin.median - bin.stdev).clamp(0.0, f32::MAX) as u32
        }) {
            let og_mdata = &itv.metadata;
            log::debug!(
                "Filtered out {:?}: {ctg}:{}-{} on contig end with coverage ({}) below bin median {:?}",
                og_mdata.typ,
                itv.first,
                itv.last,
                og_mdata.cov,
                &bin_stats[&og_mdata.bin]
            );
            *itv = Interval::new(
                itv.first,
                itv.last,
                CallInfo {
                    typ: MisassemblyType::Null,
                    cov: og_mdata.cov,
                    bin: og_mdata.bin,
                    zscore: og_mdata.zscore,
                    af: og_mdata.af,
                },
            );
            idx_end -= 1
        }
    }
}

pub(crate) fn merge_misassemblies(
    df_itvs: DataFrame,
    bin_stats: HashMap<u32, BinStats>,
    ctg: &str,
    fasta: Option<impl AsRef<Path> + Debug>,
    ignore_itvs: Option<&COITree<String, usize>>,
    cfg: Config,
) -> eyre::Result<LazyFrame> {
    let bp_merge = cfg.general.bp_merge.try_into()?;
    let cfg_min_size = cfg.minimum_size.unwrap_or_default();

    let itvs_all: Vec<(u64, u64, u32, &str, u32, f32, f32)> = multizip((
        df_itvs.column("st")?.u64()?.iter().flatten(),
        df_itvs.column("end")?.u64()?.iter().flatten(),
        df_itvs.column("cov")?.u32()?.iter().flatten(),
        df_itvs.column("status")?.str()?.iter().flatten(),
        df_itvs.column("bin")?.u32()?.iter().flatten(),
        // This is really strange behavior, the last f32 is None when using f32().iter() but a valid value in the dataframe.
        // If unhandled, multizip will omit the last record.
        // I don't have a clean solution here. This might produce a false zscore at the end of the window.
        df_itvs
            .column("zscore")?
            .f32()?
            .iter()
            .map(|v| v.unwrap_or_default()),
        df_itvs
            .column("af")?
            .f32()?
            .iter()
            .map(|v| v.unwrap_or_default()),
    ))
    .collect();
    // crate::io::write_tsv(&mut df_itvs.clone(), Some("test.tsv"))?;

    let (Some(all_st), Some(all_end)) = (
        itvs_all.first().map(|itv| itv.0 as i32),
        itvs_all.last().map(|itv| itv.1 as i32),
    ) else {
        bail!("No intervals for {ctg}. Something is wrong.");
    };

    let df_misasm_itvs = df_itvs
        .clone()
        .lazy()
        .filter(col("status").neq(lit("correct")))
        .collect()?;

    // TODO: Rewrite merging function to merge over three intervals
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
        |a, b| a.metadata.0 == b.metadata.0,
        |itv_1, itv_2| (itv_1.metadata.0, (itv_1.metadata.1 + itv_2.metadata.1) / 2),
        |itv| itv,
    );
    let final_misasm_itvs: COITree<(MisassemblyType, u32), usize> =
        COITree::new(itvs_misasm.iter());
    let thr_minimum_sizes: HashMap<MisassemblyType, u64> = (&cfg_min_size).try_into()?;

    let mut fasta_reader = if let Some(fasta) = fasta {
        log::info!("Reading indexed {fasta:?} for {ctg} to classify misassemblies by repeat.");
        Some(FastaHandle::new(fasta)?)
    } else {
        None
    };

    // Convert good intervals overlapping misassembly types.
    // Detect repeats.
    // Remove ignored intervals.
    let mut reclassified_itvs_all: Vec<Interval<CallInfo>> = Vec::with_capacity(itvs_all.len());
    for (st, end, cov, status, bin, zscore, af) in itvs_all {
        let st = st.try_into()?;
        let end = end.try_into()?;
        let len = (end - st) as f32;
        let mut largest_ovl: Option<MisassemblyType> = None;
        let mtype = MisassemblyType::from_str(status)?;
        final_misasm_itvs.query(st, end, |ovl_itv| {
            let ovl_prop = overlap_length(st, end, ovl_itv.first, ovl_itv.last) as f32 / len;
            // Needs majority.
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
            .unwrap_or(mtype);

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
            let end = end
                .saturating_add(cfg_rpt.bp_extend.try_into()?)
                .try_into()?;
            let record = reader.fetch(ctg, st.try_into()?, end)?;
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
                reclassified_itvs_all.push(Interval::new(
                    itv.first,
                    itv.last,
                    CallInfo {
                        typ: status,
                        cov,
                        bin,
                        zscore: OrderedFloat(zscore),
                        af: OrderedFloat(af),
                    },
                ));
            }
            continue;
        }

        // Otherwise, add misassembly.
        reclassified_itvs_all.push(Interval::new(
            st,
            end,
            CallInfo {
                typ: status,
                cov,
                bin,
                zscore: OrderedFloat(zscore),
                af: OrderedFloat(af),
            },
        ));
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
    let mut minmax_reclassified_itvs_all: Vec<Interval<CallInfo>> = vec![];
    for ((is_mergeable, bin), group_elements) in &reclassified_itvs_all
        .into_iter()
        .chunk_by(|a| (a.metadata.typ.is_mergeable(), a.metadata.bin))
    {
        if is_mergeable {
            let (mut agg_st, mut agg_end, mut mean_cov, mut max_zscore, mut max_af) = (
                i32::MAX,
                0,
                0,
                OrderedFloat(f32::MIN),
                OrderedFloat(f32::MIN),
            );
            let mut agg_status = MisassemblyType::Null;
            let mut num_elems = 0;
            // Get min max of region.
            for (st, end, status, cov, zscore, af) in group_elements
                .map(|itv| {
                    (
                        itv.first,
                        itv.last,
                        itv.metadata.typ,
                        itv.metadata.cov,
                        itv.metadata.zscore,
                        itv.metadata.af,
                    )
                })
                .sorted_by(|a, b| a.0.cmp(&b.0))
            {
                agg_st = std::cmp::min(st, agg_st);
                agg_end = std::cmp::max(end, agg_end);
                agg_status = std::cmp::max(agg_status, status);
                max_zscore = std::cmp::max(max_zscore, zscore);
                max_af = std::cmp::max(max_af, af);
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
            // Can possibly have no bin if entire region is misassembled.
            if let Some(bin_stats) = bin_stats.get(&bin) {
                if agg_status == MisassemblyType::Collapse
                    && bin_stats
                        .itree_above_median
                        .query_count(agg_st, agg_end)
                        .ge(&1)
                {
                    log::debug!("Filtered out {agg_status:?}: {ctg}:{agg_st}-{agg_end} above median coverage at bin transition ({bin_stats:?})");
                    agg_status = MisassemblyType::Null;
                }
            }

            minmax_reclassified_itvs_all.push(Interval::new(
                agg_st,
                agg_end,
                CallInfo {
                    typ: agg_status,
                    cov: mean_cov,
                    bin,
                    zscore: max_zscore,
                    af: max_af,
                },
            ));
        } else {
            minmax_reclassified_itvs_all.extend(group_elements.into_iter().map(|itv| {
                Interval::new(
                    itv.first,
                    itv.last,
                    CallInfo {
                        typ: itv.metadata.typ,
                        cov: itv.metadata.cov,
                        bin: itv.metadata.bin,
                        zscore: itv.metadata.zscore,
                        af: itv.metadata.af,
                    },
                )
            }));
        }
    }

    let fn_finalizer = |a: Interval<CallInfo>| -> Interval<CallInfo> {
        let mut status = a.metadata.typ;
        // Remove misassemblies less than threshold size.
        let min_size = thr_minimum_sizes[&status];
        let length = (a.last - a.first) as u64;
        if length < min_size {
            status = MisassemblyType::Null;
        }
        Interval::new(
            a.first,
            a.last,
            CallInfo {
                typ: status,
                cov: a.metadata.cov,
                bin: a.metadata.bin,
                zscore: a.metadata.zscore,
                af: a.metadata.af,
            },
        )
    };
    // Remove intervals not within minimum sizes after merging.
    // Then, remerge intervals.
    let minmax_reclassified_itvs_all = merge_intervals(
        minmax_reclassified_itvs_all.into_iter(),
        1,
        |a, b| a.metadata.typ == b.metadata.typ,
        |a, b| CallInfo {
            typ: a.metadata.typ,
            cov: (a.metadata.cov + b.metadata.cov) / 2,
            bin: a.metadata.bin,
            zscore: a.metadata.zscore.max(b.metadata.zscore),
            af: a.metadata.af.max(b.metadata.af),
        },
        fn_finalizer,
    );

    let minmax_reclassified_itvs_all: Vec<Row> = minmax_reclassified_itvs_all
        .into_iter()
        .map(|itv| {
            Row::new(vec![
                AnyValue::Int32(itv.first),
                AnyValue::Int32(itv.last),
                AnyValue::String(if itv.metadata.typ == MisassemblyType::Null {
                    "correct"
                } else {
                    itv.metadata.typ.into()
                }),
                AnyValue::UInt32(itv.metadata.cov),
                AnyValue::Float32(*itv.metadata.zscore),
                AnyValue::Float32(*itv.metadata.af),
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
            ("zscore".into(), DataType::Float32),
            ("af".into(), DataType::Float32),
        ]),
    )?;

    // Reduce final interval groups to min/max.
    Ok(df_itvs_all
        .lazy()
        .with_column(col("status").rle_id().alias("group"))
        .group_by(["group"])
        .agg([
            // Positions in noodles are 1-based ([start, end]) so need to shift.
            // https://github.com/zaeleus/noodles/discussions/207
            col("st").min() - lit(1),
            col("end").max() - lit(1),
            col("cov").median().cast(DataType::UInt32),
            col("status").first(),
            col("zscore").max(),
            col("af").max(),
        ])
        .with_column(
            when(col("status").eq(lit("correct")))
                .then(lit(0.0))
                .otherwise(col("zscore"))
                .alias("zscore"),
        )
        .sort(["st"], Default::default())
        .select([
            col("st"),
            col("end"),
            col("status"),
            col("cov"),
            col("zscore"),
            col("af"),
        ]))
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
    let thr_false_dup = (cfg.cov.ratio_false_dup * median_cov as f32).floor();
    let thr_collapse = (cfg.cov.ratio_collapse * median_cov as f32).floor();

    let lf_pileup = lf_pileup
        .with_columns([
            (col("insertion").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                .alias("insertion_ratio"),
            (col("deletion").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                .alias("deletion_ratio"),
            (col("softclip").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                .alias("softclip_ratio"),
            (col("mismatch").cast(DataType::Float32) / col("cov").cast(DataType::Float32))
                .alias("mismatch_ratio"),
        ])
        .with_column(
            // indels
            // Region with insertion or deletion peak and has high indel ratio
            //
            // deletion
            when(
                col("deletion_ratio")
                    .gt_eq(lit(cfg.indel.ratio_deletion))
                    .and(col("deletion_peak").eq(lit("high")))
                    .and(col("deletion").gt_eq(col("insertion"))),
            )
            .then(lit("deletion"))
            // insertion
            .when(
                col("insertion_ratio")
                    .gt_eq(lit(cfg.indel.ratio_insertion))
                    .and(col("insertion_peak").eq(lit("high")))
                    .and(col("insertion").gt_eq(col("deletion"))),
            )
            .then(lit("insertion"))
            // softclip
            .when(
                col("softclip_ratio")
                    .gt_eq(lit(cfg.softclip.ratio_softclip))
                    .and(col("softclip_peak").eq(lit("high"))),
            )
            .then(lit("softclip"))
            .otherwise(lit("correct"))
            .alias("status"),
        )
        .with_column(
            // collapse
            // Regions with at double the coverage and increase in mismatches/indels.
            when(
                col("cov_peak")
                    .eq(lit("high"))
                    .and(col("cov").gt_eq(lit(thr_collapse)))
                    .and(
                        col("mismatch_peak")
                            .eq(lit("high"))
                            .or(col("insertion_peak").eq(lit("high")))
                            .or(col("deletion_peak").eq(lit("high"))),
                    ),
            )
            .then(lit("collapse"))
            // misjoin
            // Regions with zero coverage.
            .when(col("cov").eq(lit(0)))
            .then(lit("misjoin"))
            // false_dup
            // Region with half of the expected coverage and a maximum mapq of zero due to multiple primary mappings.
            // Either a duplicated contig, duplicated region, or an SV (large insertion of repetive region).
            .when(
                col("cov")
                    .lt_eq(lit(thr_false_dup))
                    .and(col("mapq_max").eq(lit(0))),
            )
            .then(lit("false_dup"))
            .otherwise(col("status"))
            .alias("status"),
        )
        .with_column(
            // mismatch
            // Region that mismatches the assembly and has non-zero coverage.
            when(col("cov").eq(col("mismatch")).and(col("cov").neq(lit(0))))
                .then(lit("mismatch"))
                // het_mismap
                // Regions with high mismatch peak and het ratio.
                .when(
                    col("mismatch_ratio")
                        .gt_eq(lit(cfg.mismatch.ratio_het))
                        .and(col("mismatch_peak").eq(lit("high"))),
                )
                .then(lit("het_mismap"))
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
            // Only use correct regions for bin stats.
            .filter(col("status").eq(lit("correct")))
            .with_columns([
                col("cov")
                    .rolling_median(rolling_opts.clone())
                    .alias("cov_median"),
                col("cov").rolling_std(rolling_opts).alias("cov_stdev"),
            ])
            .select([col("bin"), col("cov_median"), col("cov_stdev")])
            .collect()?;

        // Don't use polars ChunkedArray::first() as unchecked and will segfault if empty despite type signature being Option<T>
        // https://docs.rs/polars-core/0.50.0/src/polars_core/chunked_array/mod.rs.html#568
        let bin = df_bin_stats
            .column("bin")?
            .u32()?
            .iter()
            .flatten()
            .next()
            .unwrap_or_default();
        BinStats {
            num: bin,
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
    let cols = [
        col("chrom"),
        col("pos"),
        col("cov"),
        col("status"),
        col("mismatch"),
        col("mapq"),
        col("insertion"),
        col("deletion"),
        col("softclip"),
        col("bin"),
        col("bin_ident"),
        col("zscore"),
        col("af"),
    ];
    let df_pileup = lf_pileup
        .with_columns([
            lit(ctg).alias("chrom"),
            // Get the z-score.
            when(col("status").eq(lit("deletion")))
                .then(col("deletion_zscore"))
                .when(col("status").eq(lit("insertion")))
                .then(col("insertion_zscore"))
                .when(
                    col("status")
                        .eq(lit("mismatch"))
                        .or(col("status").eq(lit("het_mismap"))),
                )
                .then(col("mismatch_zscore"))
                .when(col("status").eq(lit("collapse")))
                .then(col("cov_zscore"))
                .when(col("status").eq(lit("softclip")))
                .then(col("softclip_zscore"))
                .otherwise(lit(0.0))
                .alias("zscore"),
            // Collect the allele frequencies.
            when(col("status").eq(lit("deletion")))
                .then(col("deletion_ratio"))
                .when(col("status").eq(lit("insertion")))
                .then(col("insertion_ratio"))
                .when(
                    col("status")
                        .eq(lit("mismatch"))
                        .or(col("status").eq(lit("het_mismap"))),
                )
                .then(col("mismatch_ratio"))
                .when(col("status").eq(lit("collapse")))
                .then(max_horizontal([
                    col("insertion_ratio").max(),
                    col("deletion_ratio").max(),
                    col("mismatch_ratio").max(),
                ])?)
                .when(col("status").eq(lit("softclip")))
                .then(col("softclip_ratio"))
                .otherwise(lit(0.0))
                .alias("af"),
        ])
        .select(cols)
        .collect()?;

    // Construct intervals.
    // Store [st,end,type,cov,status,bin,zscore]
    let df_itvs = df_pileup
        .select(["pos", "cov", "status", "bin", "zscore", "af"])?
        .lazy()
        .with_column(
            ((col("pos") - col("pos").shift_and_fill(1, 0))
                .lt_eq(1)
                .rle_id()
                + col("status").rle_id())
            .alias("group"),
        );
    let df_itvs = df_itvs
        .group_by([col("group")])
        .agg([
            col("pos").min().alias("st"),
            col("pos").max().alias("end") + lit(1),
            col("cov").mean().cast(DataType::UInt32),
            col("status").first(),
            col("bin").first(),
            col("zscore").max(),
            col("af").max(),
        ])
        .drop(Selector::ByName {
            names: Arc::new(["group".into()]),
            strict: true,
        })
        .collect()?;

    Ok((df_itvs, df_pileup, bin_stats))
}
