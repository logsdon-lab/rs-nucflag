use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    str::FromStr,
};

use coitrees::Interval;
use itertools::Itertools;
use polars::prelude::*;

use crate::{config::Config, preset::Preset};

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
/// use rs_nucflag::io::read_bed;
/// use coitrees::Interval;
///
/// let records = read_bed(
///     "core/test/standard/region.bed",
///     |name: &str, start: u64, stop: u64, other_cols: &str| Interval::new(start as i32, stop as i32, None::<&str>)
/// );
/// ```
/// BED4 record
/// ```
/// use rs_nucflag::io::read_bed;
/// use coitrees::Interval;
///
/// let records = read_bed(
///     "core/test/standard/region.bed",
///     |name: &str, start: u64, stop: u64, other_cols: &str| Interval::new(start as i32, stop as i32, Some(other_cols.to_owned()))
/// );
/// ```
pub fn read_bed<T: Clone + Debug>(
    bed: impl AsRef<Path>,
    intervals_fn: impl Fn(&str, u64, u64, &str) -> Interval<T>,
) -> Option<Vec<Interval<T>>> {
    let mut intervals = Vec::new();
    let bed_fh = File::open(bed).ok()?;
    let bed_reader = BufReader::new(bed_fh);

    for line in bed_reader.lines() {
        let Ok(line) = line else {
            log::error!("Invalid line: {line:?}");
            continue;
        };
        let (name, start, stop, other_cols) =
            if let Some((name, start, stop, other_cols)) = line.splitn(4, '\t').collect_tuple() {
                (name, start, stop, other_cols)
            } else if let Some((name, start, stop)) = line.splitn(3, '\t').collect_tuple() {
                (name, start, stop, "")
            } else {
                log::error!("Invalid line: {line}");
                continue;
            };
        let (Ok(first), Ok(last)) = (start.parse::<u64>(), stop.parse::<u64>()) else {
            log::error!("Cannot parse {start} or {stop} for {line}");
            continue;
        };

        intervals.push(intervals_fn(name, first, last, other_cols))
    }
    Some(intervals)
}

pub fn read_cfg(path: Option<impl AsRef<Path>>, preset: Option<&str>) -> eyre::Result<Config> {
    match (path, preset.map(Preset::from_str)) {
        (None, None) => Ok(Config::default()),
        (None, Some(preset)) => {
            let preset = preset?;
            Ok(Config::from(preset))
        }
        (Some(cfg_path), None) => {
            let cfg_str = std::fs::read_to_string(cfg_path)?;
            toml::from_str(&cfg_str).map_err(Into::into)
        }
        (Some(cfg_path), Some(preset)) => {
            let preset = preset?;
            let cfg_str = std::fs::read_to_string(cfg_path)?;
            let cfg: Config = toml::from_str(&cfg_str)?;
            let preset_cfg = Config::from(preset);
            Ok(cfg.merge(preset_cfg))
        }
    }
}
