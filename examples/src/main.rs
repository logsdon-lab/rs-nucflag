use std::error::Error;

use coitrees::Interval;
use rayon::{prelude::*, ThreadPoolBuilder};
use rs_nucflag::{
    classify::NucFlagResult,
    config::Config,
    io::{read_bed, read_cfg, write_tsv},
    nucflag,
};

/*
cargo run --release --manifest-path examples/Cargo.toml -- \
core/test/ending_scaffold/aln_1.bam \
core/test/ending_scaffold/aln_1.fa \
core/test/ending_scaffold/aln_1.bed \
core/nucflag.toml \
1
*/
fn cli_fasta() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 6,
        "Usage: .{} <bam> <fasta> <bed> <cfg> <threads>",
        args[0]
    );

    let bam = args.get(1).expect("No bamfile provided.");
    let fasta = args.get(2).filter(|fasta| fasta.as_str() != "none");
    let bed = args.get(3).unwrap();
    let config = args.get(4);
    let threads = args
        .get(5)
        .map(|arg| arg.parse::<usize>())
        .unwrap_or(Ok(1))?;
    let cfg: Config = read_cfg(config, None)?;

    // Set number of threads.
    ThreadPoolBuilder::new().num_threads(threads);

    let ctg_itvs: Vec<Interval<String>> = read_bed(bed, |name, st, end, _| {
        Interval::new(
            st.try_into().unwrap(),
            end.try_into().unwrap(),
            name.to_owned(),
        )
    })
    .unwrap();

    // Parallelize by contig.
    let all_regions: Vec<NucFlagResult> = ctg_itvs
        .into_par_iter()
        .map(|itv| {
            // Open the BAM file in read-only per thread.
            nucflag(bam, fasta, &itv, None, cfg.clone()).unwrap()
        })
        .collect();

    eprintln!("Done! Writing to stdout.");
    for mut res in all_regions.into_iter() {
        write_tsv(&mut res.regions, None::<&str>)?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    cli_fasta()?;
    Ok(())
}
