use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use coitrees::Interval;
use itertools::Itertools;
use polars::prelude::*;

use crate::config::Config;

/// Write TSV file to file or stdout.
pub fn write_tsv(df: &mut DataFrame, path: Option<impl AsRef<Path>>) -> eyre::Result<()> {
    let mut file: Box<dyn Write> = if let Some(path) = path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(df)?;
    Ok(())
}

#[allow(unused)]
pub fn write_itvs<'a, T: Debug + Clone + 'a>(
    itvs: impl Iterator<Item = Interval<&'a T>>,
    path: Option<impl AsRef<Path>>,
) -> eyre::Result<()> {
    let mut file: Box<dyn Write> = if let Some(path) = path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    for itv in itvs {
        writeln!(&mut file, "{}\t{}\t{:?}", itv.first, itv.last, itv.metadata)?;
    }
    Ok(())
}

/// Read a BED file and return a list of [`Interval`]s.
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
///     |name: &str, start: u64, stop: u64, other_cols: &str| Interval::new(start, stop, None)
/// )
/// ```
/// BED4 record
/// ```
/// let records = read_bed(
///     "test.bed",
///     |name: &str, start: u64, stop: u64, other_cols: &str| Interval::new(start, stop, Some(other_cols.to_owned()))
/// )
/// ```
pub fn read_bed<T: Clone + Debug>(
    bed: Option<impl AsRef<Path>>,
    intervals_fn: impl Fn(&str, u64, u64, &str) -> Interval<T>,
) -> eyre::Result<Vec<Interval<T>>> {
    let mut intervals = Vec::new();

    let Some(bed) = bed else {
        return Ok(intervals);
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
        let (first, last) = (start.parse::<u64>()?, stop.parse::<u64>()?);

        intervals.push(intervals_fn(name, first, last, other_cols))
    }
    Ok(intervals)
}

pub fn read_cfg(path: Option<impl AsRef<Path>>) -> eyre::Result<Config> {
    if let Some(cfg_path) = path {
        let cfg_str = std::fs::read_to_string(cfg_path)?;
        toml::from_str(&cfg_str).map_err(Into::into)
    } else {
        Ok(Config::default())
    }
}
