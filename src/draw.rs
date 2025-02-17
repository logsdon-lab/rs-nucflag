use coitrees::Interval;
use eyre::bail;
use plotters::{
    prelude::*,
    style::full_palette::{ORANGE, PURPLE, TEAL},
};
use polars::prelude::*;

use crate::{misassembly::Misassembly, peak::Peak};

pub fn draw_nucfreq(
    outfile: &str,
    dim: (u32, u32),
    name: &str,
    positions: &Column,
    first_data: &Column,
    second_data: &Column,
    mapq_data: &Column,
    first_peaks: Vec<Interval<Peak>>,
    second_peaks: Vec<Interval<Peak>>,
) -> eyre::Result<()> {
    let root_area = BitMapBackend::new(outfile, dim).into_drawing_area();

    root_area.fill(&WHITE)?;

    let root_area = root_area.titled(name, ("sans-serif", 34))?;
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
    let Some(mean_first_cov) = first_data.u64()?.max() else {
        bail!("No data to calculate mean for {name}.")
    };

    let x_range = min_pos..max_pos + 1;
    let y_range = 0u64..mean_first_cov as u64;

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
        first_data
            .u64()?
            .iter()
            .flatten()
            .zip(x_range.clone())
            .map(|(y, x)| (x, y)),
        &BLACK,
    ))?;
    cc.draw_series(LineSeries::new(
        second_data
            .u64()?
            .iter()
            .flatten()
            .zip(x_range.clone())
            .map(|(y, x)| (x, y)),
        &RED,
    ))?;

    cc.draw_series(LineSeries::new(
        mapq_data
            .u8()?
            .iter()
            .flatten()
            .zip(x_range.clone())
            .map(|(y, x)| (x, y as u64)),
        &PURPLE,
    ))?;

    cc.draw_series(first_peaks.iter().map(|pk| {
        Rectangle::new(
            [
                (pk.first.try_into().unwrap(), y_range.end),
                (pk.last.try_into().unwrap(), y_range.start),
            ],
            ShapeStyle {
                color: (if pk.metadata == Peak::Low {
                    ORANGE
                } else {
                    GREEN
                })
                .to_rgba()
                .mix(0.5),
                filled: true,
                stroke_width: 1,
            },
        )
    }))?;

    cc.draw_series(second_peaks.iter().map(|pk| {
        Rectangle::new(
            [
                (pk.first.try_into().unwrap(), y_range.end),
                (pk.last.try_into().unwrap(), y_range.start),
            ],
            ShapeStyle {
                color: TEAL.to_rgba().mix(0.5),
                filled: true,
                stroke_width: 1,
            },
        )
    }))?;

    root_area.present()?;
    Ok(())
}
