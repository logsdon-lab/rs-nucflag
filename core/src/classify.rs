use core::str;
use std::{collections::HashMap, fmt::Debug, fs::File, io::BufReader, path::Path};

use crate::{
    config::Config,
    intervals::{merge_intervals, subtract_intervals},
    pileup::PileupInfo,
    repeats::detect_largest_repeat,
};
use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use itertools::{multizip, Itertools};
use noodles::{fasta, core::{Position, Region}};
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
    let new_itvs = subtract_intervals(curr_itv, &all_ovls);

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
    ctg: &str,
    fasta: Option<impl AsRef<Path> + Debug>,
    ignore_itvs: Option<&COITree<String, usize>>,
    cfg: Config,
) -> eyre::Result<LazyFrame> {
    let bp_merge = cfg.general.bp_merge.try_into()?;
    let merge_across_type = cfg.general.merge_across_type;
    let cfg_min_size = cfg.minimum_size.unwrap_or_default();

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

    let mut fasta_reader = if let Some(fasta) = fasta {
        log::info!("Reading indexed {fasta:?} for {ctg} to detect further classify misassemblies by repeat.");
        let mut fai = fasta.as_ref().as_os_str().to_owned();
        fai.push(".fai");
        let fai = fasta::fai::read(fai)?;
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
    let mut minimum_sizes = Vec::with_capacity(itvs_all.len());
    // TODO: Parallelize.
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

        let status = if let Some((new_status, _)) = largest_ovl {
            new_status.to_owned()
        } else {
            status
        };

        // Detect scaffold/homopolymer/simple repeat and replace type.
        let status = if let Some(reader) = fasta_reader
            .as_mut()
            .filter(|_| status == "indel" || status == "misjoin")
        {
            let region = Region::new(
                ctg,
                Position::new(st.clamp(1, i32::MAX).try_into()?).unwrap()
                    ..=Position::new(end.try_into()?).unwrap(),
            );
            let record = reader.query(&region)?;
            let seq = str::from_utf8(record.sequence().as_ref())?;
            detect_largest_repeat(seq)
                .and_then(|rpt| (rpt.prop > 0.5).then(|| rpt.repeat.to_string()))
                .unwrap_or(status)
        } else {
            status
        };

        // TODO: This might not be the best approach, but it's the easiest :)
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
                minimum_sizes.push(
                    thr_minimum_sizes
                        .get(status.as_str())
                        .cloned()
                        .unwrap_or(u64::MIN),
                );
                statuses.push(status.clone());
            }
            continue;
        }

        // Otherwise, add misassembly.
        sts.push(st);
        ends.push(end);
        covs.push(cov);
        minimum_sizes.push(
            thr_minimum_sizes
                .get(status.as_str())
                .cloned()
                .unwrap_or(u64::MIN),
        );
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
