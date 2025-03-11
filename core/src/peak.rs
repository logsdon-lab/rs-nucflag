use eyre::bail;
use polars::prelude::*;

// Polars port of https://github.com/swizard0/smoothed_z_score/blob/master/src/lib.rs
// Remove influence parameter.
pub fn find_peaks(
    df_pileup: DataFrame,
    n_zscore_low: f32,
    n_zscore_high: f32,
    min_perc: Option<f32>,
) -> eyre::Result<LazyFrame> {
    assert_eq!(
        df_pileup.get_column_names().len(),
        2,
        "Only two columns expected (pos, data)."
    );

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
            PlHashMap::from_iter([(colname.as_str(), DataType::Float32)]),
            true,
        )
        .select([col("pos"), col(&colname)]);

    let median_col = format!("{colname}_median");
    let stdev_col = format!("{colname}_stdev");
    let peak_col = format!("{colname}_peak");
    let zscore_col = format!("{colname}_zscore");

    let lf_pileup = lf_pileup
        .with_column(col(&colname).median().alias(&median_col))
        .with_column(col(&colname).median().sqrt().alias(&stdev_col))
        // Calculate zscore.
        .with_column(((col(&colname) - col(&median_col)) / col(&stdev_col)).alias(&zscore_col))
        .with_column(
            when(col(&zscore_col).gt(lit(n_zscore_high)))
                .then(lit("high"))
                .when(col(&zscore_col).lt(lit(-n_zscore_low)))
                .then(lit("low"))
                .otherwise(lit("null"))
                .alias(&peak_col),
        );

    // Filter by percentile if added.
    let lf_pileup = if let Some(min_perc) = min_perc {
        lf_pileup
            .with_column(
                col(&colname)
                    .rank(RankOptions::default(), None)
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

    // Go back to u64
    Ok(lf_pileup.cast(
        PlHashMap::from_iter([(colname.as_str(), DataType::UInt64)]),
        true,
    ))
}
