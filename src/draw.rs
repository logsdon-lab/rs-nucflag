use std::{path::Path, str::FromStr};

use eyre::bail;
use itertools::{multizip, Itertools};
use plotters::{prelude::*, style::full_palette::PURPLE};
use polars::prelude::*;

use crate::misassembly::MisassemblyType;

pub fn draw_nucfreq(
    outfile: impl AsRef<Path>,
    dim: (u32, u32),
    title: &str,
    df_cov: &DataFrame,
    df_itvs: &DataFrame,
) -> eyre::Result<()> {
    let [positions, first, second, _, _, mapq, _] = df_cov.get_columns() else {
        bail!("Invalid number of columns. Developer error.")
    };
    let root_area = BitMapBackend::new(&outfile, dim).into_drawing_area();

    root_area.fill(&WHITE)?;

    let root_area = root_area.titled(title, ("sans-serif", 34))?;
    // https://users.rust-lang.org/t/mini-rfc-min-max/19995/4
    let Some((min_pos, max_pos)) = positions
        .u64()?
        .iter()
        .flatten()
        .fold(None, |m: Option<(u64, u64)>, x| {
            m.map_or(Some((x, x)), |(m1, m2)| Some((m1.min(x), m2.max(x))))
        })
    else {
        bail!("No min or max position.")
    };

    // Find average in first_data.
    // TODO: Make configurable.
    let Some(mean_first_cov) = first.u64()?.max() else {
        bail!("No data to calculate mean for {title}.")
    };

    let x_range = min_pos..max_pos + 1;
    let y_range = 0..mean_first_cov;

    let mut cc = ChartBuilder::on(&root_area)
        .margin(5)
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(x_range.clone(), y_range.clone())?;

    cc.configure_mesh()
        .x_desc("Position")
        .y_desc("Coverage")
        .disable_mesh()
        .draw()?;

    cc.draw_series(LineSeries::new(
        first
            .u64()?
            .iter()
            .flatten()
            .zip(x_range.clone())
            .map(|(y, x)| (x, y)),
        &BLACK,
    ))?;
    cc.draw_series(LineSeries::new(
        second
            .u64()?
            .iter()
            .flatten()
            .zip(x_range.clone())
            .map(|(y, x)| (x, y)),
        &RED,
    ))?;

    // TODO: Make colorscale.
    // cc.draw_series(LineSeries::new(
    //     mapq.u8()?
    //         .iter()
    //         .flatten()
    //         .zip(x_range.clone())
    //         .map(|(y, x)| (x, y as u64)),
    //     &PURPLE,
    // ))?;

    let (sts, ends, statuses): (&Column, &Column, &Column) = df_itvs
        .columns(["st", "end", "status"])?
        .into_iter()
        .collect_tuple()
        .unwrap();

    cc.draw_series(
        multizip((
            sts.u64()?.iter().flatten(),
            ends.u64()?.iter().flatten(),
            statuses.str()?.iter().flatten(),
        ))
        .map(|(st, end, status)| {
            Rectangle::new(
                [(st, y_range.end), (end, y_range.start)],
                ShapeStyle {
                    color: RGBAColor::from(MisassemblyType::from_str(status).unwrap()),
                    filled: true,
                    stroke_width: 1,
                },
            )
        }),
    )?;

    root_area.present()?;
    Ok(())
}
