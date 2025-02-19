use std::fmt::Debug;

use eyre::bail;
use polars::prelude::*;

use crate::Interval;

pub type Position<T> = (u64, T);

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
    let lf_pileup = df_pileup
        .lazy()
        // Cast to i64 as need negative values to detect both peaks/valleys
        .cast(
            PlHashMap::from_iter([(colname.as_str(), DataType::Int64)]),
            true,
        )
        .select([col("pos"), col(&colname)]);

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

    let mean_col = format!("{colname}_mean");
    let stdev_col = format!("{colname}_stdev");

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
        .with_column(zscore_peak_call!(colname, mean_col, stdev_col, n_zscore))
        // Go back to i64
        .cast(
            PlHashMap::from_iter([(colname.as_str(), DataType::UInt64)]),
            true,
        );
    Ok(lf_pileup)
}

/// Merge positions `(u64, T)` into [`Interval`]s.
///
/// # Args
/// * `positions`: Positions iterator with variable metadata type. Must be sorted.
/// * `fn_cmp_eq`: Function to determine equality and whether to group positions.
/// * `fn_build_itv`: Function to construct interval.
/// * `fn_bp_merge`: Distance in base pairs to merge interval.
/// * `thr_bp_peak`: Minimum size of interval.
///
pub fn merge_positions<T: Debug + Clone>(
    positions: impl Iterator<Item = Position<T>>,
    fn_cmp_eq: impl Fn(&Position<T>, &Position<T>) -> bool,
    fn_build_itv: impl Fn(&[Position<T>]) -> eyre::Result<Interval<T>>,
    bp_merge: Option<u64>,
    min_bp_itv: Option<u64>,
) -> eyre::Result<Vec<Interval<T>>> {
    let mut final_itvs: Vec<Interval<T>> = vec![];
    let mut curr_itv: Vec<Position<T>> = vec![];
    let bp_merge = bp_merge.unwrap_or(1);
    let min_bp_itv = min_bp_itv.unwrap_or(1);

    for pos in positions.into_iter() {
        if let Some(prev_pos) = curr_itv.last() {
            // If separate peak type, create new.
            if fn_cmp_eq(prev_pos, &pos) {
                let dst_diff = pos.0 - prev_pos.0;
                // Greater than merge distance, add current peak elements. Then store new peak.
                if dst_diff > bp_merge {
                    final_itvs.push(fn_build_itv(&curr_itv)?);
                    curr_itv.clear();
                }
            } else {
                final_itvs.push(fn_build_itv(&curr_itv)?);
                curr_itv.clear()
            }
        }
        curr_itv.push(pos);
    }
    if !curr_itv.is_empty() {
        final_itvs.push(fn_build_itv(&curr_itv)?);
        curr_itv.clear();
    }
    final_itvs.retain(|itv| itv.length() >= min_bp_itv);
    Ok(final_itvs)
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use crate::Interval;

    use super::merge_positions;

    fn assert_itvs_equal<T>(itvs_1: impl Iterator<Item = T>, itvs_2: impl Iterator<Item = T>)
    where
        T: PartialEq + Debug,
    {
        itvs_1
            .into_iter()
            .zip(itvs_2.into_iter())
            .for_each(|(i1, i2)| assert_eq!(i1, i2));
    }

    fn build_interval<'a>(
        curr_peak: &[(u64, (&'a str, u64))],
    ) -> eyre::Result<Interval<(&'a str, u64)>> {
        let curr_peak_len: u64 = curr_peak.len().try_into()?;
        // Find min and max position. Also calculate arithmetic mean.
        let (pk_first, pk_last, pk_value_sum) =
            curr_peak.iter().fold((u64::MAX, 0u64, 0u64), |acc, x| {
                let (st, end, val_sum) = acc;
                let (nxt_pos, (_, nxt_val)) = x;
                (
                    std::cmp::min(st, *nxt_pos),
                    std::cmp::max(end, *nxt_pos),
                    val_sum + *nxt_val,
                )
            });
        let curr_peak_type = curr_peak
            .first()
            .map(|(_, pk_mdata)| pk_mdata.0)
            .unwrap_or("null");
        Ok(Interval::new(
            pk_first.try_into()?,
            pk_last.try_into()?,
            (curr_peak_type, pk_value_sum / curr_peak_len),
        ))
    }

    #[test]
    fn test_merge_peaks() {
        // (pos, (type, value))
        let peaks: Vec<(u64, (&str, u64))> = vec![
            (1, ("null", 1)),
            (2, ("null", 1)),
            (3, ("null", 1)),
            (4, ("null", 1)),
            (7, ("null", 2)),
            (8, ("null", 4)),
            (1222, ("null", 1)),
        ];
        let itvs = merge_positions(
            peaks.into_iter(),
            |m1, m2| m1 == m2,
            build_interval,
            Some(1),
            Some(0),
        )
        .unwrap();

        assert_itvs_equal(
            itvs.into_iter().map(|itv| (itv.st, itv.end, itv.metadata)),
            [
                (1, 4, ("null", 1)),
                (7, 8, ("null", 3)),
                (1222, 1222, ("null", 1)),
            ]
            .into_iter(),
        )
    }

    #[test]
    fn test_merge_sep_peaks() {
        // (pos, (type, value))
        let peaks: Vec<(u64, (&str, u64))> = vec![
            (1, ("null", 1)),
            (2, ("null", 1)),
            (3, ("high", 1)),
            (4, ("low", 1)),
            (7, ("null", 2)),
            (8, ("null", 4)),
            (1222, ("null", 1)),
        ];
        let itvs = merge_positions(
            peaks.into_iter(),
            |m1, m2| m1 == m2,
            build_interval,
            Some(1),
            Some(0),
        )
        .unwrap();

        assert_itvs_equal(
            itvs.into_iter().map(|itv| (itv.st, itv.end, itv.metadata)),
            [
                (1, 2, ("null", 1)),
                (3, 3, ("high", 1)),
                (4, 4, ("low", 1)),
                (7, 8, ("null", 3)),
                (1222, 1222, ("null", 1)),
            ]
            .into_iter(),
        )
    }
}
