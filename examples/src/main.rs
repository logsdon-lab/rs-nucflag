use std::error::Error;

use coitrees::Interval;
use nucflag::{
    classify::{nucflag, NucFlagResult},
    config::Config,
    io::{read_bed, read_cfg, write_tsv},
};
use rayon::{prelude::*, ThreadPoolBuilder};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 5,
        "Usage: .{} <bam> <bed> <cfg> <threads>",
        args[0]
    );

    let bam = args.get(1).expect("No bamfile provided.");
    let bed = args.get(2);
    let config = args.get(3);
    let threads = args
        .get(4)
        .map(|arg| arg.parse::<usize>())
        .unwrap_or(Ok(1))?;
    let cfg: Config = read_cfg(config)?;

    // Set number of threads.
    ThreadPoolBuilder::new().num_threads(threads);

    let ctg_itvs: Vec<Interval<String>> = read_bed(bed, |name, st, end, _| {
        Interval::new(
            st.try_into().unwrap(),
            end.try_into().unwrap(),
            name.to_owned(),
        )
    })?;

    // Parallelize by contig.
    let all_regions: Vec<NucFlagResult> = ctg_itvs
        .into_par_iter()
        .map(|itv| {
            // Open the BAM file in read-only per thread.
            nucflag(bam, &itv, cfg.clone(), None).unwrap()
        })
        .collect();

    eprintln!("Done! Writing to stdout.");
    for mut res in all_regions.into_iter() {
        write_tsv(&mut res.regions, None::<&str>)?;
    }

    Ok(())
}
