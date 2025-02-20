use eyre::bail;
use polars::prelude::*;

macro_rules! zscore_peak_call {
    ($colname:ident, $colname_mean:ident, $colname_stdev:ident, $colname_peak:ident, $n_zscore:expr) => {
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
        .otherwise(lit("null"))
        .alias(&$colname_peak)
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

    let mean_col = format!("{colname}_mean");
    let stdev_col = format!("{colname}_stdev");
    let peak_col = format!("{colname}_peak");

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
        .with_column(zscore_peak_call!(
            colname, mean_col, stdev_col, peak_col, n_zscore
        ))
        // Go back to i64
        .cast(
            PlHashMap::from_iter([(colname.as_str(), DataType::UInt64)]),
            true,
        );

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
            .with_column(
                when((col("rank") / col("rank").max()).gt(lit(min_perc)))
                    .then(col(&peak_col))
                    .otherwise(lit("null")),
            )
    } else {
        lf_pileup
    };

    Ok(lf_pileup)
}

#[cfg(test)]
mod test {}
