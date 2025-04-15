use core::str;
use std::{collections::HashMap, path::Path};

use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;
use noodles::{
    core::{Position, Region},
    fasta,
};
use ordered_float::OrderedFloat;
use polars::prelude::*;
use rs_moddotplot::{compute_local_seq_self_identity, compute_seq_self_identity};

pub fn split_pileup(
    mut df: DataFrame,
    fasta: impl AsRef<Path>,
    itv: &Interval<String>,
    window_size: usize,
    thr_dt_ident: f32,
) -> eyre::Result<DataFrame> {
    let ctg = itv.metadata.clone();
    let (st, end): (i32, i32) = (itv.first, itv.last);

    let mut reader_fasta = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;

    let position = Position::new(itv.first.try_into()?).unwrap()
        ..=Position::new(itv.last.try_into()?).unwrap();
    let region = Region::new(itv.metadata.clone(), position);
    let seq = reader_fasta.query(&region)?;

    log::info!("Calculating self-identity for {ctg}:{st}-{end} to bin region.");
    let itv_idents: COITree<u64, usize> = {
        let bed_ident = compute_seq_self_identity(
            str::from_utf8(seq.sequence().as_ref())?,
            &itv.metadata,
            Some(rs_moddotplot::SelfIdentConfig {
                window_size,
                ..Default::default()
            }),
        );
        let bed_local_ident = compute_local_seq_self_identity(
            &bed_ident,
            Some(rs_moddotplot::LocalSelfIdentConfig {
                window_size,
                ..Default::default()
            }),
        );
        // Group ident values.
        // We use ordered floats with treat NaNs differently. We should never have a NaN.
        let mut all_idents = bed_local_ident
            .iter()
            .map(|row| OrderedFloat(row.avg_perc_id_by_events))
            .sorted_by(|a, b| a.cmp(b))
            .dedup_by(|a, b| a == b)
            .peekable();

        let mut groups: HashMap<OrderedFloat<f32>, u64> = HashMap::new();
        let mut group_n: u64 = 0;
        let thr_dt_ident = OrderedFloat(thr_dt_ident);
        while let Some(ident) = all_idents.next() {
            let Some(next_ident) = all_idents.peek() else {
                group_n += 1;
                groups.insert(ident, group_n);
                break;
            };
            groups.insert(ident, group_n);

            if *next_ident - ident > thr_dt_ident {
                group_n += 1;
            }
        }
        log::info!("Split {ctg}:{st}-{end} into {group_n} region(s).");

        COITree::new(
            &bed_local_ident
            .into_iter()
            .map(|row|
                // Adjust coordinates.
                Interval::new(row.start as i32 + st, row.end as i32 + st, groups[&OrderedFloat(row.avg_perc_id_by_events)])
            ).collect::<Vec<Interval<u64>>>()
        )
    };

    // Add groups to pileup.
    // N's will cause offset so need to detect overlaps.
    let ident_groups: Vec<u64> = df
        .column("pos")?
        .cast(&DataType::Int32)?
        .i32()?
        .into_iter()
        .flatten()
        .map(|p| {
            let mut ident = None;
            itv_idents.query(p, p + 1, |itv| ident = Some(itv.metadata));
            ident.unwrap_or_default()
        })
        .collect();

    df.with_column(Column::new("bin".into(), ident_groups))?;

    Ok(df)
}
