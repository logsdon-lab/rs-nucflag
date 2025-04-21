use core::str;
use std::{path::Path, str::FromStr};

use crate::{
    binning::group_pileup_by_ani,
    classify::{classify_peaks, merge_misassemblies, merge_pileup_info, NucFlagResult},
    config::Config,
    misassembly::MisassemblyType,
    peak::find_peaks,
    pileup::AlignmentFile,
};
use coitrees::{COITree, Interval};
use polars::prelude::*;
use rayon::prelude::*;

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
    let lf_cov_peaks = lf_cov_peaks.join(
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
    let lf_pileup = lf_cov_peaks
        .join(
            df_pileup.select(["pos", "mapq"])?.lazy(),
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        .with_column(
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
/// * `ignore_itvs`: Intervals to ignore. NOTE: Interval's metadata is not checked.
/// * `cfg`: Peak-calling configuration. See [`Preset`] for configuration based on sequencing data type.
///
/// # Returns
/// * [`NucFlagResult`]
pub fn nucflag(
    aln: impl AsRef<Path>,
    fasta: Option<impl AsRef<Path> + Clone>,
    itv: &Interval<String>,
    ignore_itvs: Option<&COITree<String, usize>>,
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
        ignore_itvs,
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
