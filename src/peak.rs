use std::{convert::Infallible, str::FromStr};

use coitrees::{GenericInterval, Interval};
use eyre::bail;
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

macro_rules! zscore_peak_call {
    ($colname:ident, $colname_mean:ident, $colname_stdev:ident, $n_zscore:expr) => {
        when(
            // (|x - mean| / stdev) > stdev * n_stdev
            ((col(&$colname) - col(&$colname_mean)) / col(&$colname_stdev))
                .abs()
                .gt(lit($n_zscore)),
        )
        .then(
            // Then classify is peak or valley
            when(col(&$colname).gt(col(&$colname_mean)))
                .then(lit("high"))
                .otherwise(lit("low")),
        )
        .otherwise(lit("none"))
        .alias(format!("{}_peak", &$colname))
    };
}

// Polars port of https://github.com/swizard0/smoothed_z_score/blob/master/src/lib.rs
// Remove influence parameter.
pub fn find_peaks(
    df_pileup: DataFrame,
    window_size: usize,
    n_zscore: f32,
    min_perc: Option<f32>,
) -> eyre::Result<LazyFrame> {
    assert_eq!(
        df_pileup.get_column_names().len(),
        2,
        "Only two columns expected (pos, data)."
    );

    let window_opts = RollingOptionsFixedWindow {
        window_size,
        center: true,
        ..Default::default()
    };

    let Some(colname) = df_pileup
        .get_column_names_str()
        .into_iter()
        .find(|c| (*c).ne("pos"))
        .map(|c| c.to_owned())
    else {
        bail!("No colname.")
    };
    let df_pileup = df_pileup
        .lazy()
        .cast(
            PlHashMap::from_iter([(colname.as_str(), DataType::Int64)]),
            true,
        )
        .collect()?;
    let lf_pileup = df_pileup.select(["pos", &colname])?.lazy();
    // Filter by percentile if added.
    let lf_pileup = if let Some(min_perc) = min_perc {
        lf_pileup
            .with_column(
                col(&colname)
                    .rank(RankOptions::default(), None)
                    .cast(DataType::Float32)
                    .alias("rank"),
            )
            // Filter positions under threshold percentile
            // Removes noise and take only most prominent peaks.
            .filter((col("rank") / col("rank").max()).gt(lit(min_perc)))
    } else {
        lf_pileup
    };

    // // TODO: Might need to look at mapq for HSAT and false dupes.
    let mean_col = format!("{}_mean", colname);
    let stdev_col = format!("{}_stdev", colname);

    let lf_pileup = lf_pileup
        .with_column(
            col(&colname)
                .rolling_mean(window_opts.clone())
                .fill_null(col(&colname).mean())
                .alias(&mean_col),
        )
        .with_column(
            col(&colname)
                .rolling_std(window_opts.clone())
                .fill_null(col(&colname).std(1))
                .alias(&stdev_col),
        )
        .with_column(zscore_peak_call!(colname, mean_col, stdev_col, n_zscore));

    Ok(lf_pileup)
}

pub fn merge_peaks(
    peaks: impl Iterator<Item = (u64, Peak)>,
    bp_merge: Option<u64>,
    thr_bp_peak: Option<u64>,
) -> eyre::Result<Vec<Interval<Peak>>> {
    let mut final_peaks: Vec<Interval<Peak>> = vec![];
    let mut curr_peak: Vec<(u64, Peak)> = vec![];
    let bp_merge = bp_merge.unwrap_or(1);
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
