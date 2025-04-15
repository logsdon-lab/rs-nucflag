use coitrees::Interval;
use core::str;
use nucflag::{
    classify::nucflag,
    io::{read_bed, read_cfg},
    pileup::AlignmentFile,
};
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};

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

/// Classify a missassembly.
#[pyfunction]
#[pyo3(signature = (aln, fasta = None, bed = None, threads = 1, cfg = None))]
fn run_nucflag(
    aln: &str,
    fasta: Option<&str>,
    bed: Option<&str>,
    threads: usize,
    cfg: Option<&str>,
) -> PyResult<Vec<PyNucFlagResult>> {
    let cfg = read_cfg(cfg).map_err(|err| PyValueError::new_err(err.to_string()))?;

    if cfg.general.verbose {
        simple_logger::init_with_level(log::Level::Debug).expect("Cannot initialize logger.");
    }

    // Set rayon threadpool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    let ctg_itvs: Vec<Interval<String>> = if bed.is_none() {
        // If no intervals, apply to whole genome based on header.
        let mut aln =
            AlignmentFile::new(aln, fasta).map_err(|err| PyValueError::new_err(err.to_string()))?;
        let header = aln
            .header()
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        let window = cfg.general.bp_wg_window;
        header
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
            .collect()
    } else {
        read_bed(bed, |name, st, end, _| {
            Interval::new(
                st.try_into().unwrap(),
                end.try_into().unwrap(),
                name.to_owned(),
            )
        })
        .map_err(|err| PyValueError::new_err(err.to_string()))?
    };

    // Parallelize by contig.
    Ok(ctg_itvs
        .into_par_iter()
        .flat_map(|itv| {
            // Open the BAM file in read-only per thread.
            let res = nucflag(aln, fasta, &itv, cfg.clone());
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
