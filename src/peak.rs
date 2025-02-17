use std::{convert::Infallible, str::FromStr};

use coitrees::{GenericInterval, Interval};
use polars::prelude::*;

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
            (col($colname) - col(concat!($colname, "_median")))
                .abs()
                .gt(col(concat!($colname, "_stdev")) * lit($n_stdevs)),
        )
        .then(
            when(col($colname).gt(col(concat!($colname, "_median"))))
                .then(lit("high"))
                .otherwise(lit("low")),
        )
        .otherwise(lit("none"))
        .alias(concat!($colname, "_peak"))
    };
}

pub fn find_peaks(df_pileup: DataFrame, n_stdevs: i64) -> eyre::Result<DataFrame> {
    let window_opts = RollingOptionsFixedWindow {
        window_size: 5000,
        ..Default::default()
    };
    let lf_pileup = df_pileup.lazy().cast(
        PlHashMap::from_iter([
            ("first", DataType::Int64),
            ("second", DataType::Int64),
            ("mapq", DataType::Int8),
        ]),
        true,
    );
    let lf_pileup = lf_pileup
        .with_columns([
            col("first")
                .rolling_median(window_opts.clone())
                .alias("first_median"),
            col("second")
                .rolling_median(window_opts.clone())
                .alias("second_median"),
            col("mapq")
                .rolling_median(window_opts.clone())
                .alias("mapq_median"),
            col("first")
                .rolling_std(window_opts.clone())
                .alias("first_stdev"),
            col("second")
                .rolling_std(window_opts.clone())
                .alias("second_stdev"),
            col("mapq")
                .rolling_std(window_opts.clone())
                .alias("mapq_stdev"),
        ])
        .with_columns([
            expr_peaks!("first", n_stdevs),
            expr_peaks!("second", n_stdevs),
            expr_peaks!("mapq", n_stdevs),
        ])
        .select([
            col("pos"),
            col("first"),
            col("second"),
            col("mapq"),
            col("sec"),
            col("sup"),
            col("first_peak"),
            col("second_peak"),
            col("mapq_peak"),
        ])
        .cast(
            PlHashMap::from_iter([
                ("first", DataType::UInt64),
                ("second", DataType::UInt64),
                ("mapq", DataType::UInt8),
            ]),
            true,
        );

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

    fn build_interval(curr_peak: &Vec<(u64, Peak)>) -> eyre::Result<Interval<Peak>> {
        let (curr_peak_first, curr_peak_last) =
            (curr_peak.first().unwrap(), curr_peak.last().unwrap());
        Ok(Interval::new(
            curr_peak_first.0.try_into()?,
            curr_peak_last.0.try_into()?,
            curr_peak_first.1,
        ))
    };

    for pk in peaks {
        if let Some(prev_pk) = curr_peak.last().filter(|prev_pk| prev_pk.1 == pk.1) {
            let dst_diff = pk.0 - prev_pk.0;
            // Greater than merge distance, add current peak elements. Then store new peak.
            if dst_diff > bp_merge {
                final_peaks.push(build_interval(&curr_peak)?);
                curr_peak.clear();
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
