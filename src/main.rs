use clap::Parser;
use classify::classify_misassemblies;
use cli::Cli;
use coitrees::Interval;
use io::{read_bed, read_cfg, write_tsv};
use noodles::bam::{self};
use polars::prelude::*;
use rayon::{prelude::*, ThreadPoolBuilder};

mod classify;
mod cli;
mod config;
mod draw;
mod intervals;
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
        Interval::new(st, end, name.to_owned())
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
    let all_bed_misassemblies: Vec<LazyFrame> = ctg_itvs
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

    let mut all_misassemblies = concat(
        all_bed_misassemblies,
        UnionArgs {
            parallel: true,
            ..Default::default()
        },
    )?
    .collect()?;

    // Write all misassemblies to file or stdout
    write_tsv(&mut all_misassemblies, outfile)?;

    log::info!("Done!");
    Ok(())
}
