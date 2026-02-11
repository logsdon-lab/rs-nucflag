use coitrees::{COITree, Interval, IntervalTree};
use core::str;
use pyo3::{exceptions::PyValueError, prelude::*};
use rs_nucflag::{io::read_bed, pileup::AlignmentFile};
use std::collections::HashMap;

pub(crate) fn get_whole_genome_intervals(
    aln: &str,
    window: usize,
) -> Result<Vec<Interval<String>>, PyErr> {
    // If no intervals, apply to whole genome based on header.
    let mut aln = AlignmentFile::new(aln).map_err(|err| PyValueError::new_err(err.to_string()))?;
    aln.aligned_intervals_windows(window)
        .map_err(|err| PyValueError::new_err(err.to_string()))
}

pub(crate) fn get_aln_intervals(
    aln: &str,
    bed: Option<&str>,
    bp_wg_window: usize,
) -> Result<Vec<Interval<String>>, PyErr> {
    let all_intervals = get_whole_genome_intervals(aln, bp_wg_window)
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    if let Some(bed) = bed {
        let all_intervals_map: HashMap<String, Interval<String>> = all_intervals
            .into_iter()
            .map(|itv| (itv.metadata.clone(), itv))
            .collect();
        let valid_bed_intervals = read_bed(bed, |name, st, end, _| {
            Interval::new(st as i32, end as i32, name.to_owned())
        })
        .map(|itvs| {
            itvs.into_iter()
                .filter_map(|itv| {
                    if let Some(itv_seq_bounds) = all_intervals_map.get(&itv.metadata) {
                        let new_first = itv.first.clamp(0, itv_seq_bounds.last);
                        let new_last = itv.last.clamp(0, itv_seq_bounds.last);
                        Some(Interval::new(new_first, new_last, itv.metadata))
                    } else {
                        log::error!("Skipping invalid interval: {itv:?}");
                        None
                    }
                })
                .collect()
        })
        .ok_or_else(|| PyValueError::new_err(format!("Unable to read intervals from {bed}")))?;
        Ok(valid_bed_intervals)
    } else {
        Ok(all_intervals)
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
