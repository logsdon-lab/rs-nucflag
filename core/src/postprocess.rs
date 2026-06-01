use std::{collections::HashMap, fmt::Debug, path::Path, str::FromStr};

use crate::{
    binning::BinStats,
    config::Config,
    intervals::{merge_intervals, overlap_length},
    io::FastaHandle,
    misassembly::MisassemblyType,
    repeats::{detect_largest_repeat, Repeat},
};
use coitrees::{COITree, Interval, IntervalTree};
use eyre::bail;
use itertools::{multizip, Itertools};
use ordered_float::OrderedFloat;
use polars::{frame::row::Row, prelude::*};

#[derive(Debug, Clone)]
pub(crate) struct CallInfo {
    typ: MisassemblyType,
    cov: u32,
    bin: u32,
    zscore: OrderedFloat<f32>,
    af: OrderedFloat<f32>,
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

pub(crate) fn postprocess_misassemblies(
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
    );
    let final_misasm_itvs: COITree<(MisassemblyType, u32), usize> =
        COITree::new(itvs_misasm.iter());
    let thr_minimum_sizes: HashMap<MisassemblyType, u64> = (&cfg_min_size).try_into()?;
    // crate::io::write_itvs(itvs_misasm.iter().cloned(), Some("test_before.tsv"))?;

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
        let mut status = if let (Some(reader), Some(cfg_rpt)) = (
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
        if ignore_itvs.is_some_and(|itree| itree.query_count(st, end) > 0) {
            status = MisassemblyType::Null
        }

        // Handle case where should not add insertion/deletion between repeat types
        if matches!(
            status,
            MisassemblyType::Insertion | MisassemblyType::Deletion
        ) && zscore == 0.0
        {
            status = mtype
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
    reclassified_itvs_all.sort_by_key(|a| a.first);
    // crate::io::write_itvs(reclassified_itvs_all.iter().cloned(), Some("test_after.tsv"))?;

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

    // Remove intervals not within minimum sizes after merging.
    // Then, remerge intervals.
    let mut minmax_reclassified_itvs_all = merge_intervals(
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
    );

    // Remove misassemblies less than threshold size.
    for itv in minmax_reclassified_itvs_all.iter_mut() {
        let min_size = thr_minimum_sizes[&itv.metadata.typ];
        let length = (itv.last - itv.first) as u64;
        if length < min_size {
            itv.metadata.typ = MisassemblyType::Null;
        }
    }

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
