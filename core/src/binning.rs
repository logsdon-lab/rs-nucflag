use core::str;
use std::{collections::VecDeque, fmt::Debug, path::Path};

use coitrees::{COITree, Interval, IntervalTree};
use polars::prelude::*;
use rs_moddotplot::{compute_group_seq_self_identity, compute_seq_self_identity};

use crate::{config::GroupByANIConfig, io::FastaHandle};

pub struct BinStats {
    pub num: u32,
    pub median: f32,
    pub stdev: f32,
    pub itree_above_median: COITree<(), usize>,
}

impl Debug for BinStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BinStats")
            .field("num", &self.num)
            .field("median", &self.median)
            .field("stdev", &self.stdev)
            .finish()
    }
}

macro_rules! get_median_cov {
    ($lf:expr, $start:expr, $end:expr) => {
        $lf.filter(col("pos").gt_eq($start).and(col("pos").lt_eq($end)))
            .select([col("cov")])
            .collect()?
            .column("cov")?
            .median_reduce()?
            .value()
            .try_extract::<u32>()
            .unwrap_or_default()
    };
}

pub fn group_pileup_by_ani(
    mut df: DataFrame,
    fasta: impl AsRef<Path> + Debug,
    itv: &Interval<String>,
    cfg: &GroupByANIConfig,
) -> eyre::Result<DataFrame> {
    let ctg = itv.metadata.clone();
    let (st, end): (i32, i32) = (itv.first, itv.last);
    let st_print = st - (itv.first == 1) as i32;
    let window_size = cfg.window_size;
    let min_grp_size = cfg.min_grp_size;
    let min_ident = cfg.min_ident;

    let mut reader_fasta = FastaHandle::new(fasta)?;
    let seq = reader_fasta.fetch(&ctg, st as usize, end as usize)?;

    let itv_idents: COITree<(u32, f32), usize> = {
        log::info!(
            "Calculating self-identity for {ctg}:{}-{end} to bin region.",
            st_print
        );
        let bed_ident = compute_seq_self_identity(
            str::from_utf8(seq.sequence().as_ref())?,
            &itv.metadata,
            Some(rs_moddotplot::SelfIdentConfig {
                window_size,
                ..Default::default()
            }),
        );
        log::info!("Grouping repetitive intervals in {ctg}:{}-{end}.", st_print);
        let bed_group_ident = compute_group_seq_self_identity(&bed_ident);

        let mut itvs: VecDeque<Interval<(u32, f32)>> = bed_group_ident
            .into_iter()
            .filter(|r| (r.end - r.start > min_grp_size) && r.avg_perc_id_by_events > min_ident)
            .enumerate()
            // 0 is unbinned.
            .map(|(i, r)| {
                Interval::new(
                    r.start as i32 + st,
                    r.end as i32 + st,
                    ((i + 1) as u32, r.avg_perc_id_by_events),
                )
            })
            .collect();

        let mut final_itvs: Vec<Interval<(u32, f32)>> = vec![];
        let lf_cov = df.select(["pos", "cov"])?.lazy();

        while let Some(mut itv) = itvs.pop_front() {
            let Some(mut next_itv) = itvs.pop_front() else {
                final_itvs.push(itv);
                break;
            };
            // Ensure no small gaps between.
            let dst_between = next_itv.first - itv.last;
            if dst_between <= min_grp_size.try_into()? {
                let cov_left = get_median_cov!(lf_cov.clone(), itv.first, itv.last) as i32;
                let cov_between = get_median_cov!(lf_cov.clone(), itv.last, next_itv.first) as i32;
                let cov_right =
                    get_median_cov!(lf_cov.clone(), next_itv.first, next_itv.last) as i32;

                let cov_diff_left = cov_between.abs_diff(cov_left);
                let cov_diff_right = cov_between.abs_diff(cov_right);

                // Calculate cov difference and choose side to fill gap.
                if cov_diff_left < cov_diff_right {
                    itv.last = next_itv.first;
                } else {
                    next_itv.first = itv.last;
                }
                itv.last = next_itv.first;
            }
            final_itvs.push(itv);
            // Add to front to continue merging.
            itvs.push_front(next_itv);
        }
        COITree::new(&final_itvs)
    };
    log::info!(
        "Detected {} region(s) in {ctg}:{}-{end}.",
        itv_idents.len() + 1,
        st_print,
    );

    // Add groups to pileup.
    // N's will cause offset so need to detect overlaps.
    let (ident_groups, ident_values): (Vec<u32>, Vec<f32>) = df
        .column("pos")?
        .cast(&DataType::Int32)?
        .i32()?
        .into_iter()
        .flatten()
        .map(|p| {
            let mut group = None;
            itv_idents.query(p, p + 1, |itv| group = Some(itv.metadata));
            // If not in group, assign to other bin, 0.
            group.unwrap_or_default()
        })
        .unzip();

    df.with_column(Column::new("bin".into(), ident_groups))?;
    df.with_column(Column::new("bin_ident".into(), ident_values))?;

    Ok(df)
}
