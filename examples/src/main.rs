use std::{
    fs::File,
    io::{BufWriter, Write},
};

use cli::Cli;
use nucflag::{
    classify::{classify_misassemblies, NucFlagResult},
    io::{read_bed, read_cfg},
};
use clap::Parser;
use coitrees::Interval;
use noodles::bam::{self};
use rayon::{prelude::*, ThreadPoolBuilder};

mod cli;

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
    let all_regions: Vec<NucFlagResult> = ctg_itvs
        .into_par_iter()
        .map(|itv| {
            // Open the BAM file in read-only per thread.
            classify_misassemblies(bamfile.clone(), &itv, cfg.clone(), None).unwrap()
        })
        .collect();

    // let mut output_fh: Box<dyn Write> =
    //     if let Some(outfile) = outfile.as_ref().and_then(|f| f.to_str()) {
    //         log::info!("Done! Writing to {outfile}.");
    //         let file = File::create(outfile)?;
    //         Box::new(BufWriter::new(file))
    //     } else {
    //         log::info!("Done! Writing to stdout.");
    //         let buffer = BufWriter::new(std::io::stdout());
    //         Box::new(buffer)
    //     };

    // for res in all_regions.into_iter() {
    // }

    Ok(())
}