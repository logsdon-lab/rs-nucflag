use std::{convert::Infallible, str::FromStr};

use coitrees::{GenericInterval, Interval};
use polars::prelude::*;

use crate::io::write_tsv;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Peak {
    Low,
    High,
    Null,
}

impl FromStr for Peak {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "low" => Ok(Peak::Low),
            "high" => Ok(Peak::High),
            _ => Ok(Peak::Null),
        }
    }
}

macro_rules! expr_peaks {
    ($colname:literal, $n_stdevs:expr) => {
        when(
            // |x - median| > stdev * n_stdev
            (col($colname) - col(concat!($colname, "_median")))
                .abs()
                .gt(col(concat!($colname, "_stdev")) * lit($n_stdevs)),
        )
        .then(
            // Then classify is peak or valley
            when(col($colname).gt(col(concat!($colname, "_median"))))
                .then(lit("high"))
                .otherwise(lit("low")),
        )
        .otherwise(lit("none"))
        .alias(concat!($colname, "_peak"))
    };
}

// Polars port of https://github.com/swizard0/smoothed_z_score/blob/master/src/lib.rs
// Remove influence parameter.
pub fn find_peaks(
    df_pileup: DataFrame,
    window_size: usize,
    n_stdevs: i64,
    thr_second_perc: f32,
) -> eyre::Result<DataFrame> {
    let window_opts = RollingOptionsFixedWindow {
        window_size,
        ..Default::default()
    };
    let df_pileup = df_pileup
        .lazy()
        .cast(
            PlHashMap::from_iter([("first", DataType::Int64), ("second", DataType::Int64)]),
            true,
        )
        .collect()?;
    let lf_first = df_pileup.select(["pos", "first"])?.lazy();
    let lf_second = df_pileup.select(["pos", "second"])?.lazy();
    // // TODO: Might need to look at mapq for HSAT and false dupes.

    // let lf_rest = df_pileup
    //     .lazy()
    //     .select([col("pos"), col("mapq"), col("sec"), col("sup")]);

    let lf_first = lf_first
        .with_columns([
            col("first")
                .rolling_median(window_opts.clone())
                .alias("first_median"),
            (col("first").rolling_sum(window_opts.clone()) / lit(window_size as u64))
                .sqrt()
                .alias("first_stdev"),
        ])
        .with_column(expr_peaks!("first", n_stdevs));

    // Take only second cols that are >= 50th percentile.
    let lf_second = lf_second
        .with_column(
            col("second")
                .rank(RankOptions::default(), None)
                .cast(DataType::Float32)
                .alias("rank"),
        )
        // Filter positions under threshold percentile
        // Removes noise and take only most prominent peaks.
        .filter((col("rank") / col("rank").max()).gt(lit(thr_second_perc)))
        .with_columns([
            col("second")
                .rolling_median(window_opts.clone())
                .alias("second_median"),
            (col("second")
                .rolling_sum(window_opts.clone())
                .floor_div(lit(window_size as u64)))
            .sqrt()
            .alias("second_stdev"),
        ])
        .with_column(expr_peaks!("second", n_stdevs));

    let lf_pileup = lf_first
        .join(
            lf_second,
            [col("pos")],
            [col("pos")],
            JoinArgs::new(JoinType::Left),
        )
        // .join(
        //     lf_rest,
        //     [col("pos")],
        //     [col("pos")],
        //     JoinArgs::new(JoinType::Left),
        // )
        .select([
            col("pos"),
            // col("first"),
            // col("second"),
            // col("mapq"),
            // col("sec"),
            // col("sup"),
            col("first_peak"),
            col("second_peak"),
        ])
        .with_columns([
            col("second_peak").fill_null(lit("none"))
        ]);

    Ok(lf_pileup.collect()?)
}

pub fn merge_peaks(
    peaks: impl Iterator<Item = (u64, Peak)>,
    bp_merge: Option<u64>,
    thr_bp_peak: Option<u64>,
) -> eyre::Result<Vec<Interval<Peak>>> {
    let mut final_peaks: Vec<Interval<Peak>> = vec![];
    let mut curr_peak: Vec<(u64, Peak)> = vec![];
    let bp_merge = bp_merge.unwrap_or(5000);
    let thr_bp_peak = thr_bp_peak.unwrap_or(1).try_into()?;

    fn build_interval(curr_peak: &[(u64, Peak)]) -> eyre::Result<Interval<Peak>> {
        let (curr_peak_first, curr_peak_last) =
            (curr_peak.first().unwrap(), curr_peak.last().unwrap());
        Ok(Interval::new(
            curr_peak_first.0.try_into()?,
            curr_peak_last.0.try_into()?,
            curr_peak_first.1,
        ))
    }

    for pk in peaks {
        if let Some(prev_pk) = curr_peak.last() {
            // If separate peak type, create new.
            if prev_pk.1 == pk.1 {
                let dst_diff = pk.0 - prev_pk.0;
                // Greater than merge distance, add current peak elements. Then store new peak.
                if dst_diff > bp_merge {
                    final_peaks.push(build_interval(&curr_peak)?);
                    curr_peak.clear();
                }
            } else {
                final_peaks.push(build_interval(&curr_peak)?);
                curr_peak.clear()
            }
        }
        curr_peak.push(pk);
    }
    if !curr_peak.is_empty() {
        final_peaks.push(build_interval(&curr_peak)?);
        curr_peak.clear();
    }
    final_peaks.retain(|pk| pk.len() > thr_bp_peak);
    Ok(final_peaks)
}

#[cfg(test)]
mod test {
    use super::{merge_peaks, Peak};

    #[test]
    fn test_merge_peaks() {
        let peaks: Vec<(u64, Peak)> = vec![
            (1, Peak::Null),
            (2, Peak::Null),
            (3, Peak::Null),
            (4, Peak::Null),
            (7, Peak::Null),
            (8, Peak::Null),
            (1222, Peak::Null),
        ];
        let itvs = merge_peaks(peaks.into_iter(), Some(1), Some(0)).unwrap();

        assert_eq!(
            itvs.into_iter()
                .map(|itv| (itv.first, itv.last, itv.metadata))
                .collect::<Vec<_>>(),
            vec![
                (1, 4, Peak::Null),
                (7, 8, Peak::Null),
                (1222, 1222, Peak::Null),
            ]
        )
    }

    #[test]
    fn test_merge_sep_peaks() {
        let peaks: Vec<(u64, Peak)> = vec![
            (1, Peak::Null),
            (2, Peak::Null),
            (3, Peak::High),
            (4, Peak::Low),
            (7, Peak::Null),
            (8, Peak::Null),
            (1222, Peak::Null),
        ];
        let itvs = merge_peaks(peaks.into_iter(), Some(1), Some(0)).unwrap();

        assert_eq!(
            itvs.into_iter()
                .map(|itv| (itv.first, itv.last, itv.metadata))
                .collect::<Vec<_>>(),
            vec![
                (1, 2, Peak::Null),
                (3, 3, Peak::High),
                (4, 4, Peak::Low),
                (7, 8, Peak::Null),
                (1222, 1222, Peak::Null),
            ]
        )
    }
}
