use coitrees::{COITree, Interval, IntervalTree};
use core::str;
use nucflag::{
    io::{read_bed, read_cfg},
    nucflag,
    pileup::AlignmentFile,
};
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};
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
    /// Coverage regions.
    #[pyo3(get)]
    pub cov: Option<PyDataFrame>,
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

/// Classify a missassembly.
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

    let ctg_itvs: Vec<Interval<String>> = if bed.is_none() {
        // If no intervals, apply to whole genome based on header.
        get_whole_genome_intervals(aln, fasta, cfg.general.bp_wg_window)?
    } else {
        read_bed(bed, |name, st, end, _| {
            Interval::new(st as i32, end as i32, name.to_owned())
        })
        .map_err(|err| PyValueError::new_err(err.to_string()))?
    };

    let ignore_itvs: HashMap<String, COITree<String, usize>> = {
        read_bed(ignore_bed, |name, start, stop, _| {
            Interval::new(start as i32, stop as i32, name.to_owned())
        })
        .map_err(|err| PyValueError::new_err(err.to_string()))?
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
        .collect()
    };

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
                    cov: res.cov.map(PyDataFrame),
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
    Ok(())
}
