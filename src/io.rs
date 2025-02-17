use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;
use polars::prelude::*;

pub type RegionIntervals<T> = HashMap<String, Vec<Interval<T>>>;
pub type RegionIntervalTrees<T> = HashMap<String, COITree<T, usize>>;

pub fn write_tsv(df: &mut DataFrame, path: impl AsRef<Path>) -> eyre::Result<()> {
    let mut file = File::create(path)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(df)?;
    Ok(())
}

/// Read an input bedfile and convert it to a [`COITree`].
///
/// # Arguments
/// * `bed`: Bedfile path.
/// * `intervals_fn`: Function applied to `(start, stop, other_cols)` to convert into an [`Interval`].
///
/// # Examples
/// BED3 record.
/// ```
/// let records = read_bed(
///     "test.bed",
///     |name: &str, start: i32, stop: i32, other_cols: &str| Interval::new(start, stop, None)
/// )
/// ```
/// BED4 record
/// ```
/// let records = read_bed(
///     "test.bed",
///     |name: &str, start: i32, stop: i32, other_cols: &str| Interval::new(start, stop, Some(other_cols.to_owned()))
/// )
/// ```
pub fn read_bed<T: Clone>(
    bed: Option<impl AsRef<Path>>,
    intervals_fn: impl Fn(&str, i32, i32, &str) -> Interval<T>,
) -> eyre::Result<Option<RegionIntervalTrees<T>>> {
    let mut intervals: RegionIntervals<T> = HashMap::new();
    let mut trees: RegionIntervalTrees<T> = HashMap::new();

    let Some(bed) = bed else {
        return Ok(None);
    };
    let bed_fh = File::open(bed)?;
    let bed_reader = BufReader::new(bed_fh);

    for line in bed_reader.lines() {
        let line = line?;
        let (name, start, stop, other_cols) =
            if let Some((name, start, stop, other_cols)) = line.splitn(4, '\t').collect_tuple() {
                (name, start, stop, other_cols)
            } else if let Some((name, start, stop)) = line.splitn(3, '\t').collect_tuple() {
                (name, start, stop, "")
            } else {
                log::error!("Invalid line: {line}");
                continue;
            };
        let (first, last) = (start.parse::<i32>()?, stop.parse::<i32>()?);

        intervals
            .entry(name.to_owned())
            .and_modify(|intervals| intervals.push(intervals_fn(name, first, last, other_cols)))
            .or_insert_with(|| vec![intervals_fn(name, first, last, other_cols)]);
    }
    for (roi, intervals) in intervals.into_iter() {
        trees.entry(roi).or_insert(COITree::new(&intervals));
    }
    Ok(Some(trees))
}
