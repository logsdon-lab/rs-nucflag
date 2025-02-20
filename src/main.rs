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
use itertools::{multizip, Itertools};
use misassembly::MisassemblyType;
use noodles::bam::{self};
use plotters::style::RGBAColor;

use polars::prelude::*;
use rayon::{prelude::*, ThreadPoolBuilder};
use utils::Interval;

mod classify;
mod cli;
mod config;
mod draw;
mod io;
mod misassembly;
mod peak;
mod pileup;
mod utils;

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
        Interval::new(st, end, name.to_owned()).unwrap()
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
                Interval::new(1, length, ctg_name.clone()).unwrap()
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
            let plot_fname = cli.plot_dir.as_ref().map(|plot_dir| {
                let mut plot_fname = plot_dir.join(ctg);
                plot_fname.set_extension("png");
                plot_fname
            });
            // Open the BAM file in read-only per thread.
            classify_misassemblies(bamfile.clone(), itv, cov_fname, plot_fname, cfg.clone())
                .unwrap()
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
            covs.u64()?.iter().flatten(),
            statuses.str()?.iter().flatten(),
        )) {
            let item_rgb = RGBAColor::from(MisassemblyType::from_str(status)?);
            // Write BED9
            writeln!(
                output_fh,
                "{ctg}\t{st}\t{end}\t{status:?}\t{cov}\t+\t{st}\t{end}\t{},{},{}",
                item_rgb.0, item_rgb.1, item_rgb.2
            )?;
        }
    }

    Ok(())
}
