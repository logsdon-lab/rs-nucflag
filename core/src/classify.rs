use std::fmt::Debug;

use crate::{binning::BinStats, config::Config};
use coitrees::{COITree, Interval, IntervalTree};
use polars::lazy::dsl::max_horizontal;
use polars::prelude::*;

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
                // het_or_mismap
                // Regions with high mismatch peak and exceeds min ratio but below ratio het.
                // Difficult to determine from alignments. Could be cell culture artifacts or mismapping from collapsed/misjoined sequence.
                .when(
                    col("mismatch_ratio").gt(lit(cfg.mismatch.ratio_het)).and(
                        col("mismatch_ratio")
                            .lt_eq(lit(cfg.mismatch.ratio_low_quality))
                            .and(col("mismatch_peak").eq(lit("high"))),
                    ),
                )
                .then(lit("het_or_mismap"))
                // low_quality
                // Regions with high mismatch peak and exceeds low quality ratio.
                .when(
                    col("mismatch_ratio")
                        .gt(lit(cfg.mismatch.ratio_low_quality))
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
                        .or(col("status").eq(lit("low_quality")))
                        .or(col("status").eq(lit("het_or_mismap"))),
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
                        .or(col("status").eq(lit("low_quality")))
                        .or(col("status").eq(lit("het_or_mismap"))),
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
