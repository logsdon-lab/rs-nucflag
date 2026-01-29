use coitrees::{COITree, Interval};
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};
use rs_nucflag::{io::read_cfg, nucflag};
use simple_logger::SimpleLogger;
use std::collections::HashMap;

mod utils;

use crate::utils::{get_aln_intervals, get_ignored_intervals};

/// NucFlag results.
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
    /// Pileup of region with columns:
    /// * `chrom`
    ///     * Chromosome name
    /// * `pos`
    ///     * Position
    /// * `cov`
    ///     * Coverage
    /// * `status`
    ///     * Status
    /// * `mismatch`
    ///     * Number of mismatches
    /// * `mapq`
    ///     * MAPQ
    /// * `insertion`
    ///     * Number of insertions
    /// * `deletion`
    ///     * Number of deletions
    /// * `softclip`
    ///     * Number of softclipped bases
    /// * `bin`
    ///     * Region number.
    ///     * 0 if non-repetitive and some number if repetitive
    /// * `bin_ident`
    ///     * Self-sequence identity of region.
    ///     * 0.0 if non-repetitive
    #[pyo3(get)]
    pub pileup: PyDataFrame,
    /// Regions and their status in BED9 format:
    /// * `#chrom`
    ///     * Chromosome name
    /// * `chromStart`
    ///     * Chromosome start
    /// * `chromEnd`
    ///     * Chromosome end
    /// * `name`
    ///     * Name of status.
    /// * `score`
    ///     * Coverage of region.
    /// * `strand`
    ///     * None or `.`
    /// * `thickStart`
    ///     * Same as `chromEnd`
    /// * `thickEnd`
    ///     * Same as `chromStart`
    /// * `itemRgb`
    ///     * Color of status
    /// * `zscore`
    ///     * [MAD adjusted z-score.](https://www.ibm.com/docs/en/cognos-analytics/12.0.x?topic=terms-modified-z-score)
    /// * `af`
    ///     * Allele frequency.
    ///     * `signal / coverage` where deletions are not counted toward `coverage`
    #[pyo3(get)]
    pub regions: PyDataFrame,
}

/// Get interval regions from an alignment file or bed file.
///
/// # Args
/// * `aln`
///     * Alignment as BAM or CRAM file.
/// * `bed`
///     * BED file with coordinates.
///     * If not provided, splits all regions listed in `aln` header (`@SQ`) into non-overlapping windows.
///     * Invalid intervals are ignored.
/// * `window`
///     * Window size in base pairs of non-overlapping windows if `bed` not provided.
///
/// # Returns
/// * List of intervals as tuples with the format, start, end, and chromosome name.
///
/// # Example
/// ```python
/// from py_nucflag import get_regions
///
/// # Get all possible aligned regions
/// regions = get_regions("aln.bam", window = 10_000_000)
/// # Get regions in bed file and verify that exist from alignment header.
/// regions_bed = get_regions("aln.bam", bed="regions.bed")
/// ```
#[pyfunction]
#[pyo3(signature = (aln, bed = None, window = 10_000_000))]
fn get_regions(aln: &str, bed: Option<&str>, window: usize) -> PyResult<Vec<(i32, i32, String)>> {
    // Init with default level.
    _ = SimpleLogger::new().init();

    Ok(get_aln_intervals(aln, bed, window)?
        .into_iter()
        .map(|itv| (itv.first, itv.last, itv.metadata))
        .collect())
}

/// Classify a missassembly for one interval. Identical to `run_nucflag` but only for one interval.
///
/// # Example
/// ```python
/// from py_nucflag import run_nucflag_itv
///
/// # Run on alignments to a single interval.
/// result = run_nucflag_itv(
///     "sample.bam",
///     itv=(0, 1_000_000), "chr1"),
///     fasta="sample.fa.gz"
/// )
/// print(result.regions)
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

    _ = SimpleLogger::new().with_level(cfg.general.log_level).init();

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

/// Return `NucFlag` config as TOML string from preset.
///
/// # Args
/// * `preset`
///     * `NucFlag` preset.
///     * Either `ont_r9`, `ont_r10`, or `hifi`
/// * `config`
///     * Path to configfile.
///     * If provided alongside `preset`, `preset` fields given priority.
///
/// # Returns
/// * Configuration as a TOML string.
///
/// # Example
/// ```python
/// from py_nucflag import get_config_from_preset
///
/// # Get preset parameters
/// config = get_config_from_preset(preset="hifi")
///
/// # Read in existing config and merge with preset parameters.
/// config_w_preset = get_config_from_preset(preset="ont_r10", cfg="config.toml")
#[pyfunction]
#[pyo3(signature = (preset = None, cfg = None))]
fn get_config_from_preset(preset: Option<&str>, cfg: Option<&str>) -> PyResult<String> {
    let cfg = read_cfg(cfg, preset).map_err(|err| PyValueError::new_err(err.to_string()))?;
    toml::to_string(&cfg).map_err(|err| PyValueError::new_err(err.to_string()))
}

/// Classify missassemblies from a whole-genome alignment.
///
/// # Args
/// * `aln`
///     * Alignment file as BAM or CRAM file. Requires `cs` tag if CRAM.
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
///
/// # Example
/// ```python
/// from py_nucflag import run_nucflag
///
/// # Run on all alignments.
/// # Providing a bed file is highly recommended as pileup information is stored.
/// results = run_nucflag(
///     "sample.bam",
///     bed="regions.bed",
///     fasta="sample.fa.gz"
/// )
/// # Iterate through each region's misassembly calls.
/// for result in results:
///     # The misassembly calls as BED9
///     print(result.regions)
///     # The raw pileup
///     print(result.pileup)
/// ```
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

    _ = SimpleLogger::new().with_level(cfg.general.log_level).init();

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
    m.add_function(wrap_pyfunction!(get_config_from_preset, m)?)?;
    Ok(())
}
