use coitrees::Interval;
use noodles::bam;
use nucflag::{
    classify::classify_misassemblies,
    io::{read_bed, read_cfg},
};
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use rayon::{prelude::*, ThreadPoolBuilder};

#[pyclass]
pub struct PyNucFlagResult {
    pub ctg: String,
    pub st: i32,
    pub end: i32,
    pub cov: PyDataFrame,
    pub regions: PyDataFrame,
}

/// Classify a missassembly.
#[pyfunction]
#[pyo3(signature = (bamfile, bedfile, threads, cfg))]
fn run_nucflag(
    bamfile: &str,
    bedfile: &str,
    threads: usize,
    cfg: &str,
) -> PyResult<Vec<PyNucFlagResult>> {
    let cfg = read_cfg(Some(cfg)).map_err(|err| PyValueError::new_err(err.to_string()))?;

    // Set rayon threadpool
    ThreadPoolBuilder::new().num_threads(threads);

    let bed = read_bed(Some(bedfile), |name, st, end, _| {
        Interval::new(
            st.try_into().unwrap(),
            end.try_into().unwrap(),
            name.to_owned(),
        )
    })
    .map_err(|err| PyValueError::new_err(err.to_string()))?;

    let ctg_itvs: Vec<Interval<String>> = if bed.is_empty() {
        // If no intervals, apply to whole genome based on header.
        let mut bamfile = bam::io::indexed_reader::Builder::default().build_from_path(bamfile)?;
        let header = bamfile.read_header()?;

        header
            .reference_sequences()
            .into_iter()
            .map(|(ctg, ref_seq)| {
                let ctg_name: String = ctg.clone().try_into().unwrap();
                let length = ref_seq.length().get().try_into().unwrap();
                Interval::new(1, length, ctg_name.clone())
            })
            .collect()
    } else {
        bed
    };

    // Parallelize by contig.
    Ok(ctg_itvs
        .into_par_iter()
        .map(|itv| {
            // Open the BAM file in read-only per thread.
            let res = classify_misassemblies(bamfile, &itv, cfg.clone()).unwrap();
            PyNucFlagResult {
                ctg: itv.metadata,
                st: itv.first,
                end: itv.last,
                cov: PyDataFrame(res.cov),
                regions: PyDataFrame(res.regions),
            }
        })
        .collect())
}

/// A Python module implemented in Rust.
#[pymodule]
fn rs_nucflag(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNucFlagResult>()?;
    m.add_function(wrap_pyfunction!(run_nucflag, m)?)?;
    Ok(())
}
