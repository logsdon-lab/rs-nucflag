use std::{
    fs::File,
    io::{BufWriter, Write},
    str::FromStr,
};

use crate::{
    classify::classify_misassemblies,
    cli::Cli,
    io::{read_bed, read_cfg},
};
use clap::Parser;
use coitrees::Interval;
use itertools::{multizip, Itertools};
use misassembly::MisassemblyType;
use noodles::bam::{self};
use polars::prelude::*;
use rayon::{prelude::*, ThreadPoolBuilder};

mod classify;
mod cli;
mod config;
mod intervals;
// mod draw;
mod io;
mod misassembly;
mod peak;
mod pileup;

// https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
fn main() -> eyre::Result<()> {
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()?;

    let cli = Cli::parse();
    let cfg = read_cfg(cli.config)?;
    let bedfile = cli.bed;
    let bamfile = cli.bam;
    let outfile = cli.misassemblies;

    // Set rayon threadpool
    ThreadPoolBuilder::new().num_threads(cli.threads);

    let bed = read_bed(bedfile, |name, st, end, _| {
        Interval::new(
            st.try_into().unwrap(),
            end.try_into().unwrap(),
            name.to_owned(),
        )
    })?;

    let ctg_itvs: Vec<Interval<String>> = if bed.is_empty() {
        // If no intervals, apply to whole genome based on header.
        let mut bamfile = bam::io::indexed_reader::Builder::default().build_from_path(&bamfile)?;
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
    let all_regions: Vec<(String, DataFrame)> = ctg_itvs
        .into_par_iter()
        .map(|itv| {
            let ctg = itv.metadata.as_str();
            let cov_fname = cli.cov_dir.as_ref().map(|cov_dir| {
                let mut cov_fname = cov_dir.join(ctg);
                cov_fname.set_extension("cov");
                cov_fname
            });
            // Open the BAM file in read-only per thread.
            classify_misassemblies(bamfile.clone(), itv, cov_fname, cfg.clone()).unwrap()
        })
        .collect();

    let mut output_fh: Box<dyn Write> =
        if let Some(outfile) = outfile.as_ref().and_then(|f| f.to_str()) {
            log::info!("Done! Writing to {outfile}.");
            let file = File::create(outfile)?;
            Box::new(BufWriter::new(file))
        } else {
            log::info!("Done! Writing to stdout.");
            let buffer = BufWriter::new(std::io::stdout());
            Box::new(buffer)
        };

    for (ctg, itvs) in all_regions.into_iter() {
        let (sts, ends, covs, statuses): (&Column, &Column, &Column, &Column) = itvs
            .columns(["st", "end", "cov", "status"])?
            .into_iter()
            .collect_tuple()
            .unwrap();
        for (st, end, cov, status) in multizip((
            sts.u64()?.iter().flatten(),
            ends.u64()?.iter().flatten(),
            covs.f64()?.iter().flatten(),
            statuses.str()?.iter().flatten(),
        )) {
            let item_rgb = MisassemblyType::from_str(status)?.item_rgb();
            // Write BED9
            writeln!(
                output_fh,
                "{ctg}\t{st}\t{end}\t{status:?}\t{cov}\t+\t{st}\t{end}\t{item_rgb}",
            )?;
        }
    }

    Ok(())
}
