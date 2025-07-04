use eyre::bail;
use polars::prelude::*;

pub fn find_peaks(
    df_pileup: DataFrame,
    n_zscore_low: f32,
    n_zscore_high: f32,
    drop_zeroes: bool,
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

    // Remove zeroes.
    let lf_pileup = if drop_zeroes {
        lf_pileup.filter(col(&colname).neq(lit(0)))
    } else {
        lf_pileup
    };

    let median_col = format!("{colname}_median");
    let stdev_col = format!("{colname}_stdev");
    let peak_col = format!("{colname}_peak");
    let zscore_col = format!("{colname}_zscore");

    let lf_pileup = lf_pileup
        // Calculate median
        .with_column(col(&colname).median().alias(&median_col))
        // Calculate adjusted z-score from MAD
        // https://www.statisticshowto.com/median-absolute-deviation/
        // https://www.ibm.com/docs/en/cognos-analytics/11.1.x?topic=terms-modified-z-score
        // https://www.statology.org/modified-z-score/
        .with_column(((col(&colname) - col(&median_col)).abs().median()).alias(&stdev_col))
        .with_column(
            // If MAD is zero, use mean absolute deviation and scale.
            when(col(&stdev_col).eq(lit(0.0)))
                .then((col(&colname) - col(&colname).mean()).abs().mean() * lit(1.253314))
                .otherwise(col(&stdev_col) * lit(1.486))
                .alias(&stdev_col),
        )
        .with_column(((col(&colname) - col(&median_col)) / col(&stdev_col)).alias(&zscore_col))
        .with_column(
            when(col(&zscore_col).gt(lit(n_zscore_high)))
                .then(lit("high"))
                .when(col(&zscore_col).lt(lit(-n_zscore_low)))
                .then(lit("low"))
                .otherwise(lit("null"))
                .alias(&peak_col),
        );

    // Go back to u64
    Ok(lf_pileup.cast(
        PlHashMap::from_iter([(colname.as_str(), DataType::UInt64)]),
        true,
    ))
}

#[cfg(test)]
mod test {
    use polars::prelude::*;

    use super::find_peaks;

    #[test]
    fn test_find_peaks() {
        let df = df!(
            "pos" => [0, 1, 2, 3, 4],
            "first" => [30, 60, 30, 0, 30],
        )
        .unwrap();

        let df_peaks = find_peaks(df, 1.5, 1.5, false).unwrap().collect().unwrap();
        let peaks = df_peaks.column("first_peak").unwrap();
        assert_eq!(
            vec!["null", "high", "null", "low", "null"],
            peaks.str().unwrap().iter().flatten().collect::<Vec<&str>>()
        );
    }
}
