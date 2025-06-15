use coitrees::{COITree, Interval, IntervalTree};
use core::str;
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};
use rs_nucflag::{
    io::{read_bed, read_cfg},
    nucflag,
    pileup::AlignmentFile,
};
use std::collections::HashMap;

#[pyclass]
pub struct PyNucFlagResult {
    /// Name of contig.
    #[pyo3(get)]
    pub ctg: String,
    /// Start of region.
    #[pyo3(get)]
    pub st: i32,
    /// End of region.
    #[pyo3(get)]
    pub end: i32,
    /// Pileup of regions.
    #[pyo3(get)]
    pub pileup: PyDataFrame,
    /// Regions and their status.
    #[pyo3(get)]
    pub regions: PyDataFrame,
}

pub(crate) fn get_whole_genome_intervals(
    aln: &str,
    fasta: Option<&str>,
    window: usize,
) -> Result<Vec<Interval<String>>, PyErr> {
    // If no intervals, apply to whole genome based on header.
    let mut aln =
        AlignmentFile::new(aln, fasta).map_err(|err| PyValueError::new_err(err.to_string()))?;
    let header = aln
        .header()
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
    Ok(header
        .reference_sequences()
        .into_iter()
        .flat_map(|(ctg, ref_seq)| {
            let ctg_name: String = ctg.clone().try_into().unwrap();
            let length: usize = ref_seq.length().get();
            let (num, rem) = (length / window, length % window);
            let final_start = num * window;
            let final_itv = Interval::new(
                final_start as i32,
                (final_start + rem) as i32,
                ctg_name.clone(),
            );
            (1..num + 1)
                .map(move |i| {
                    Interval::new(
                        ((i - 1) * window) as i32,
                        (i * window) as i32,
                        ctg_name.clone(),
                    )
                })
                .chain([final_itv])
        })
        .collect())
}

pub(crate) fn get_aln_intervals(
    aln: &str,
    fasta: Option<&str>,
    bed: Option<&str>,
    bp_wg_window: usize,
) -> Result<Vec<Interval<String>>, PyErr> {
    if let Some(bed) = bed {
        Ok(read_bed(bed, |name, st, end, _| {
            Interval::new(st as i32, end as i32, name.to_owned())
        })
        .ok_or_else(|| PyValueError::new_err(format!("Unable to read intervals from {bed}")))?)
    } else {
        Ok(get_whole_genome_intervals(aln, fasta, bp_wg_window)
            .map_err(|err| PyValueError::new_err(err.to_string()))?)
    }
}

pub(crate) fn get_ignored_intervals(
    ignore_bed: Option<&str>,
) -> Result<HashMap<String, COITree<String, usize>>, PyErr> {
    if let Some(intervals) = ignore_bed.and_then(|bed| {
        read_bed(bed, |name, start, stop, _| {
            Interval::new(start as i32, stop as i32, name.to_owned())
        })
    }) {
        Ok(intervals
            .into_iter()
            .fold(
                HashMap::new(),
                |mut acc: HashMap<String, Vec<Interval<String>>>, x| {
                    if acc.contains_key(&x.metadata) {
                        acc.get_mut(&x.metadata).unwrap().push(x);
                    } else {
                        acc.entry(x.metadata.clone()).or_insert(vec![x]);
                    }
                    acc
                },
            )
            .into_iter()
            .map(|(rgn, itvs)| (rgn, COITree::new(&itvs)))
            .collect())
    } else {
        Ok(HashMap::default())
    }
}

/// Get interval regions from an alignment file or bed file.
#[pyfunction]
#[pyo3(signature = (aln, fasta = None, bed = None, window = 10_000_000))]
fn get_regions(
    aln: &str,
    fasta: Option<&str>,
    bed: Option<&str>,
    window: usize,
) -> PyResult<Vec<(i32, i32, String)>> {
    Ok(get_aln_intervals(aln, fasta, bed, window)?
        .into_iter()
        .map(|itv| (itv.first, itv.last, itv.metadata))
        .collect())
}

/// Classify a missassembly for one interval. Identical to `run_nucflag` but only for one interval.
#[pyfunction]
#[pyo3(signature = (aln, itv, fasta = None, ignore_bed = None, threads = 1, cfg = None, preset = None))]
fn run_nucflag_itv(
    aln: &str,
    itv: (i32, i32, String),
    fasta: Option<&str>,
    ignore_bed: Option<&str>,
    threads: usize,
    cfg: Option<&str>,
    preset: Option<&str>,
) -> PyResult<PyNucFlagResult> {
    let cfg = read_cfg(cfg, preset).map_err(|err| PyValueError::new_err(err.to_string()))?;
    let itv = Interval::new(itv.0, itv.1, itv.2);
    if cfg.general.verbose {
        simple_logger::init_with_level(log::Level::Debug).expect("Cannot initialize logger.");
    }

    // Set rayon threadpool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    let all_ignore_itvs: HashMap<String, COITree<String, usize>> =
        get_ignored_intervals(ignore_bed)?;
    let ignore_itvs = all_ignore_itvs.get(&itv.metadata);
    // Open the BAM file in read-only per thread.
    nucflag(aln, fasta, &itv, ignore_itvs, cfg.clone())
        .map(|res| PyNucFlagResult {
            ctg: itv.metadata,
            st: itv.first,
            end: itv.last,
            pileup: PyDataFrame(res.pileup),
            regions: PyDataFrame(res.regions),
        })
        .map_err(|err| PyValueError::new_err(err.to_string()))
}

#[pyfunction]
#[pyo3(signature = (preset = None, cfg = None))]
fn print_config_from_preset(preset: Option<&str>, cfg: Option<&str>) -> PyResult<()> {
    let cfg = read_cfg(cfg, preset).map_err(|err| PyValueError::new_err(err.to_string()))?;
    if cfg.general.verbose {
        simple_logger::init_with_level(log::Level::Debug).expect("Cannot initialize logger.");
    }
    log::info!("Using config:\n{cfg:#?}");
    Ok(())
}

/// Classify a missassembly from an alignment file.
///
/// # Args
/// * `aln`
///     * Alignment file as BAM or CRAM file. Requires fasta if CRAM.
/// * `bed`
///     * BED3 file with regions to evaluate.
/// * `ignore_bed`
///     * BED3 file with regions to ignore.
/// * `threads`
///     * Number of threads to spawn.
/// * `cfg`
///     * Configfile. See [`nucflag::config::Config`]
/// * `preset`
///     * Configuration for specific LR sequencing reads.
///     * Modifies `cfg` where preset specific options take priority.
///     * See [`nucflag::preset::Preset`].
///
/// # Returns
/// * A [`PyNucFlagResult`] object where:
///     * `cov` is a pileup dataframe
///     * `regions` contains all regions evaluated.
#[pyfunction]
#[pyo3(signature = (aln, fasta = None, bed = None, ignore_bed = None, threads = 1, cfg = None, preset = None))]
fn run_nucflag(
    aln: &str,
    fasta: Option<&str>,
    bed: Option<&str>,
    ignore_bed: Option<&str>,
    threads: usize,
    cfg: Option<&str>,
    preset: Option<&str>,
) -> PyResult<Vec<PyNucFlagResult>> {
    let cfg = read_cfg(cfg, preset).map_err(|err| PyValueError::new_err(err.to_string()))?;

    if cfg.general.verbose {
        simple_logger::init_with_level(log::Level::Debug).expect("Cannot initialize logger.");
    }

    log::info!("Using config:\n{cfg:#?}");

    // Set rayon threadpool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    let ctg_itvs: Vec<Interval<String>> =
        get_aln_intervals(aln, fasta, bed, cfg.general.bp_wg_window)?;

    let ignore_itvs: HashMap<String, COITree<String, usize>> = get_ignored_intervals(ignore_bed)?;

    // Parallelize by contig.
    Ok(ctg_itvs
        .into_par_iter()
        .flat_map(|itv| {
            let ignore_itvs = ignore_itvs.get(&itv.metadata);
            // Open the BAM file in read-only per thread.
            let res = nucflag(aln, fasta, &itv, ignore_itvs, cfg.clone());
            match res {
                Ok(res) => Some(PyNucFlagResult {
                    ctg: itv.metadata,
                    st: itv.first,
                    end: itv.last,
                    pileup: PyDataFrame(res.pileup),
                    regions: PyDataFrame(res.regions),
                }),
                Err(err) => {
                    log::error!("Error: {err}");
                    None
                }
            }
        })
        .collect())
}

/// NucFlag implemented in Rust.
#[pymodule]
fn py_nucflag(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNucFlagResult>()?;
    m.add_function(wrap_pyfunction!(run_nucflag, m)?)?;
    m.add_function(wrap_pyfunction!(run_nucflag_itv, m)?)?;
    m.add_function(wrap_pyfunction!(get_regions, m)?)?;
    m.add_function(wrap_pyfunction!(print_config_from_preset, m)?)?;
    Ok(())
}
