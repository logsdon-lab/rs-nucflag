use coitrees::{COITree, Interval};
use core::str;
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};
use rs_nucflag::{io::read_cfg, nucflag};
use std::collections::HashMap;
mod utils;

use crate::utils::{get_aln_intervals, get_ignored_intervals};

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

/// Get interval regions from an alignment file or bed file.
#[pyfunction]
#[pyo3(signature = (aln, bed = None, window = 10_000_000))]
fn get_regions(aln: &str, bed: Option<&str>, window: usize) -> PyResult<Vec<(i32, i32, String)>> {
    Ok(get_aln_intervals(aln, bed, window)?
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

    _ = simple_logger::init_with_level(cfg.general.log_level);

    // Set rayon threadpool
    _ = ThreadPoolBuilder::new().num_threads(threads).build_global();

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

/// Print config from preset.
#[pyfunction]
#[pyo3(signature = (preset = None, cfg = None))]
fn print_config_from_preset(preset: Option<&str>, cfg: Option<&str>) -> PyResult<()> {
    let cfg = read_cfg(cfg, preset).map_err(|err| PyValueError::new_err(err.to_string()))?;
    _ = simple_logger::init_with_level(cfg.general.log_level);
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
///     * `pileup` is a pileup dataframe
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

    _ = simple_logger::init_with_level(cfg.general.log_level);

    log::info!("Using config:\n{cfg:#?}");

    // Set rayon threadpool
    _ = ThreadPoolBuilder::new().num_threads(threads).build_global();

    let ctg_itvs: Vec<Interval<String>> = get_aln_intervals(aln, bed, cfg.general.bp_wg_window)?;
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
