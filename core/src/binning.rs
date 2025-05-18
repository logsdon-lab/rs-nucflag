use core::str;
use std::path::Path;

use coitrees::{COITree, Interval, IntervalTree};
use noodles::{
    core::{Position, Region},
    fasta,
};
use polars::prelude::*;
use rs_moddotplot::{compute_group_seq_self_identity, compute_seq_self_identity};

use crate::config::GroupByANIConfig;

#[derive(Debug, Clone)]
pub struct BinStats {
    pub num: u64,
    pub median: f32,
    pub stdev: f32,
}

pub fn group_pileup_by_ani(
    mut df: DataFrame,
    fasta: impl AsRef<Path>,
    itv: &Interval<String>,
    cfg: &GroupByANIConfig,
) -> eyre::Result<DataFrame> {
    let ctg = itv.metadata.clone();
    let (st, end): (i32, i32) = (itv.first.clamp(1, i32::MAX), itv.last);
    let window_size = cfg.window_size;
    let min_grp_size = cfg.min_grp_size;
    let min_ident = cfg.min_ident;

    let mut reader_fasta = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;
    let position = Position::new(st.try_into()?).unwrap()..=Position::new(end.try_into()?).unwrap();
    let region = Region::new(itv.metadata.clone(), position);
    let seq = reader_fasta.query(&region)?;

    let itv_idents: COITree<u64, usize> = {
        log::info!("Calculating self-identity for {ctg}:{st}-{end} to bin region.");
        let bed_ident = compute_seq_self_identity(
            str::from_utf8(seq.sequence().as_ref())?,
            &itv.metadata,
            Some(rs_moddotplot::SelfIdentConfig {
                window_size,
                ..Default::default()
            }),
        );
        log::info!("Grouping repetitive intervals in {ctg}:{st}-{end}.");
        let bed_group_ident = compute_group_seq_self_identity(&bed_ident);

        COITree::new(
            &bed_group_ident
                .into_iter()
                .filter(|r| (r.end - r.start > min_grp_size) && r.avg_perc_id_by_events > min_ident)
                .enumerate()
                // 0 is unbinned.
                .map(|(i, r)| Interval::new(r.start as i32 + st, r.end as i32 + st, (i + 1) as u64))
                .collect::<Vec<Interval<u64>>>(),
        )
    };
    log::info!(
        "Detected {} region(s) in {ctg}:{st}-{end}.",
        itv_idents.len() + 1
    );

    // Add groups to pileup.
    // N's will cause offset so need to detect overlaps.
    let ident_groups: Vec<u64> = df
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
        .collect();

    df.with_column(Column::new("bin".into(), ident_groups))?;

    Ok(df)
}
