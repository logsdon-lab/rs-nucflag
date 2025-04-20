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

use crate::config::GroupByANIConfig;

pub fn group_pileup_by_ani(
    mut df: DataFrame,
    fasta: impl AsRef<Path>,
    itv: &Interval<String>,
    cfg: &GroupByANIConfig,
) -> eyre::Result<DataFrame> {
    let ctg = itv.metadata.clone();
    let (st, end): (i32, i32) = (itv.first, itv.last);
    let window_size = cfg.window_size;
    let band_size = cfg.band_size;
    let thr_dt_ident = OrderedFloat(cfg.thr_dt_ident);
    let min_grp_size = cfg.min_grp_size;

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
                n_bins: band_size,
                ignore_bins: 1,
            }),
        );
        log::info!(
            "Grouping self-identity intervals in {ctg}:{st}-{end} by a change of {thr_dt_ident}%."
        );

        // Group ident values.
        // Sort first to ensure that sequence identiy is in ascending order. 
        // We use ordered floats with treat NaNs differently. We should never have a NaN.
        let mut all_idents = bed_local_ident
            .iter()
            .map(|row| (row.end - row.start, OrderedFloat(row.avg_perc_id_by_events)))
            .sorted_by(|a, b| a.1.cmp(&b.1))
            .peekable();

        // Group interval identities.
        let mut groups: HashMap<OrderedFloat<f32>, u64> = HashMap::new();
        let mut group_sizes: HashMap<u64, usize> = HashMap::new();
        let mut group_n: u64 = 0;
        let mut group_size: usize = 0;
        while let Some((size, ident)) = all_idents.next() {
            let Some((_, next_ident)) = all_idents.peek() else {
                groups.insert(ident, group_n);
                group_sizes.insert(group_n, group_size + size);
                break;
            };
            // Add new group
            groups.insert(ident, group_n);
            group_size += size;
            // If large jump in identity, start a new group.
            if *next_ident - ident > thr_dt_ident {
                group_sizes.insert(group_n, group_size);
                group_size = 0;
                group_n += 1;
            }
        }

        log::info!("Merging small intervals less than {min_grp_size}bp.");
        // Iterate through small groups and try to find larger group of higher identity to merge into.
        let replacement_groups: HashMap<u64, u64> = group_sizes
            .iter()
            .filter(|grp| *grp.1 < min_grp_size)
            .flat_map(|(grp, _)| {
                let mut repl_grp = *grp;
                // Find next group that is higher in identity and satifies grp_len.
                loop {
                    let grp_len = group_sizes.get(&repl_grp)?;
                    if *grp_len < min_grp_size {
                        repl_grp += 1
                    } else {
                        return Some((*grp, repl_grp));
                    }
                }
            })
            .collect();
        
        // Remove groups that don't meet grp len.
        log::info!("Split {ctg}:{st}-{end} into {group_n} region(s).");

        COITree::new(
            &bed_local_ident
                .into_iter()
                .map(|row| {
                    // Adjust coordinates.
                    let grp = groups[&OrderedFloat(row.avg_perc_id_by_events)];
                    Interval::new(
                        row.start as i32 + st,
                        row.end as i32 + st,
                        *replacement_groups.get(&grp).unwrap_or(&grp),
                    )
                })
                .collect::<Vec<Interval<u64>>>(),
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
            let mut group = None;
            itv_idents.query(p, p + 1, |itv| group = Some(itv.metadata));
            // TODO: This needs to be checked.
            group.unwrap_or_default()
        })
        .collect();

    df.with_column(Column::new("bin".into(), ident_groups))?;

    Ok(df)
}
