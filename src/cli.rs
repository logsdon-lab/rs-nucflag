use std::path::PathBuf;

use clap::Parser;

#[derive(Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    /// Input BAM file.
    #[arg(short = 'i', long)]
    pub bam: PathBuf,

    /// Input BED file.
    #[arg(short = 'b', long)]
    pub bed: Option<PathBuf>,

    /// Misassembly BED9 file.
    /// Format: contig,st,end,misassembly,cov,+,st,end,item_rgb
    #[arg(short = 'o', long)]
    pub misassemblies: Option<PathBuf>,

    /// Status BED file.
    /// Format: contig,st,end,status,{perc_good, perc_collapse, ...}
    #[arg(short = 'o', long)]
    pub summary: Option<PathBuf>,

    // Output plot directory.
    #[arg(short = 'p', long = "plot-dir")]
    pub plot_dir: Option<PathBuf>,

    /// Output coverage directory. Generates gzipped per-base coverage TSV files.
    #[arg(long)]
    pub cov_dir: Option<PathBuf>,

    /// Config file.
    #[arg(short = 'c', long)]
    pub config: Option<PathBuf>,

    /// Number of threads.
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,
}
