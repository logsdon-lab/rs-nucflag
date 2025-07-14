use core::str;
use std::{collections::HashMap, fmt::Debug, path::Path, str::FromStr};

use crate::{
    binning::{group_pileup_by_ani, BinStats},
    classify::{classify_peaks, merge_misassemblies, NucFlagResult},
    config::Config,
    misassembly::MisassemblyType,
    peak::find_peaks,
    pileup::{merge_pileup_info, AlignmentFile},
};
use coitrees::{COITree, Interval};
use polars::prelude::*;
use rayon::prelude::*;

fn nucflag_grp(
    df_pileup: DataFrame,
    cfg: &Config,
    ctg: &str,
) -> eyre::Result<(DataFrame, DataFrame, BinStats)> {
    // Calculate est coverage for region or use provided.
    let est_median_cov: u32 = df_pileup
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
        false,
    )?;
    // Call peaks in mismatch-base signal.
    let lf_mismatch_peaks = find_peaks(
        df_pileup.select(["pos", "mismatch"])?,
        cfg.mismatch.n_zscores_high,
        cfg.mismatch.n_zscores_high,
        true,
    )?;
    // Detect indel peaks.
    let lf_indel_peaks = find_peaks(
        df_pileup.select(["pos", "indel"])?,
        // Don't care about dips in indels.
        cfg.indel.n_zscores_high,
        cfg.indel.n_zscores_high,
        true,
    )?;
    // Detect softclip peaks.
    let lf_softclip_peaks = find_peaks(
        df_pileup.select(["pos", "softclip"])?,
        // Don't care about dips in indels.
        cfg.softclip.n_zscores_high,
        cfg.softclip.n_zscores_high,
        true,
    )?;
    // Detect mapq dips.
    let lf_mapq_dips = find_peaks(
        df_pileup.select(["pos", "mapq"])?,
        // Don't care about peaks in mapq.
        cfg.mapq.n_zscores_low,
        cfg.mapq.n_zscores_low,
        false,
    )?;

    let lf_pileup = lf_cov_peaks
        .join(
            lf_indel_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .with_column(col("indel_peak").fill_null(lit("null")))
        .join(
            lf_mismatch_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .with_column(col("mismatch_peak").fill_null(lit("null")))
        .join(
            lf_softclip_peaks,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .with_column(col("softclip_peak").fill_null(lit("null")))
        .join(
            lf_mapq_dips,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .join(
            df_pileup
                .select(["pos", "mapq_max", "indel", "softclip", "bin"])?
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
/// * `fasta`: Input fasta file path. Used for region binning and repeat detection.
/// * `itv`: Interval to check.
/// * `ignore_itvs`: Intervals to ignore. NOTE: Interval's metadata is not checked.
/// * `cfg`: Peak-calling configuration. See [`Preset`] for configuration based on sequencing data type.
///
/// # Returns
/// * [`NucFlagResult`]
pub fn nucflag<A, F>(
    aln: A,
    fasta: Option<F>,
    itv: &Interval<String>,
    ignore_itvs: Option<&COITree<String, usize>>,
    cfg: Config,
) -> eyre::Result<NucFlagResult>
where
    A: AsRef<Path> + Debug,
    F: AsRef<Path> + Clone + Debug,
{
    let ctg = itv.metadata.clone();
    let (st, end) = (itv.first.try_into()?, itv.last.try_into()?);

    let mut aln = AlignmentFile::new(aln)?;
    let pileup = aln.pileup(
        itv,
        cfg.indel.min_ins_size,
        cfg.indel.min_del_size,
        cfg.general.bp_min_aln_length,
    )?;

    let df_raw_pileup = merge_pileup_info(pileup.pileups, st, end, &cfg)?;
    log::info!("Detecting peaks/valleys in {ctg}:{st}-{end}.");

    let df_pileup_groups = if let (Some(fasta), Some(cfg_grp_by_ani)) = (
        fasta.clone(),
        &cfg.group_by_ani
            .as_ref()
            .filter(|cfg| cfg.window_size < (itv.last - itv.first) as usize),
    ) {
        group_pileup_by_ani(df_raw_pileup, fasta, itv, cfg_grp_by_ani)?
            .partition_by(["bin"], true)?
    } else {
        vec![df_raw_pileup
            .lazy()
            .with_column(lit(0).cast(DataType::UInt64).alias("bin"))
            .collect()?]
    };

    let (dfs_itvs, (dfs_pileup, bin_stats)): (Vec<LazyFrame>, (Vec<LazyFrame>, Vec<BinStats>)) =
        df_pileup_groups
            .into_par_iter()
            .map(|df_pileup_grp| {
                let (df_itv, df_pileup, bin_stats) =
                    nucflag_grp(df_pileup_grp, &cfg, &ctg).unwrap();
                (df_itv.lazy(), (df_pileup.lazy(), bin_stats))
            })
            .unzip();
    let df_itvs = concat(dfs_itvs, Default::default())?
        .sort(["st"], Default::default())
        .collect()?;
    let bin_stats: HashMap<u64, BinStats> = bin_stats
        .into_iter()
        .map(|bstats| (bstats.num, bstats))
        .collect();
    let df_pileup = concat(dfs_pileup, Default::default())?
        .sort(["chrom", "pos"], Default::default())
        .collect()?;

    // Then merge and filter.
    log::info!("Merging intervals in {ctg}:{st}-{end}.");
    let df_itvs_final = merge_misassemblies(df_itvs, bin_stats, &ctg, fasta, ignore_itvs, cfg)?
        .with_columns([
            lit(ctg.clone()).alias("chrom"),
            col("st").alias("thickStart"),
            col("end").alias("thickEnd"),
            lit(".").alias("strand"),
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
        pileup: df_pileup,
        regions: df_itvs_final,
    })
}
