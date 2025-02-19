use std::{
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    classify::classify_misassemblies,
    cli::Cli,
    io::{read_bed, read_cfg},
    misassembly::MisassemblyType,
};
use clap::Parser;
use noodles::bam::{self};
use plotters::style::RGBAColor;
use rayon::{prelude::*, ThreadPoolBuilder};

mod classify;
mod cli;
mod config;
mod draw;
mod io;
mod misassembly;
mod peak;
mod pileup;

#[derive(Debug, Clone)]
pub struct Interval<T: Clone + std::fmt::Debug> {
    st: u64,
    end: u64,
    metadata: T,
}

impl<T: Clone + std::fmt::Debug> Interval<T> {
    pub fn new(st: u64, end: u64, metadata: T) -> Self {
        assert!(end >= st, "Invalid interval coordinates.");
        Self { st, end, metadata }
    }

    pub fn length(&self) -> u64 {
        self.end - self.st
    }
}

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
    let all_regions: Vec<(String, Vec<Interval<(MisassemblyType, u64)>>)> = ctg_itvs
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
        for itv in itvs {
            let (status, cov) = itv.metadata;
            let (st, end) = (itv.st, itv.end);
            let item_rgb = RGBAColor::from(&status);
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
