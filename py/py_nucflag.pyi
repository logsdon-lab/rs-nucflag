import polars as pl

# flake8: noqa: PYI021

class PyNucFlagResult:
    """
    NucFlag results.
    """

    ctg: str
    """
    Name of contig.
    """
    st: int
    """
    Start of region.
    """
    end: int
    """
    End of region.
    """
    pileup: pl.DataFrame
    """
    Pileup of regions.
    """
    regions: pl.DataFrame
    """
    Regions and their status.
    """

def get_regions(
    aln: str, bed: str | None = None, window: int = 10000000
) -> list[tuple[int, int, str]]:
    """Get interval regions from an alignment file or bed file."""

def print_config_from_preset(preset: str | None = None, cfg: str | None = None) -> None:
    """Print config from preset."""

def run_nucflag(
    aln: str,
    fasta: str | None = None,
    bed: str | None = None,
    ignore_bed: str | None = None,
    threads: int = 1,
    cfg: str | None = None,
    preset: str | None = None,
) -> list[PyNucFlagResult]:
    """
    Classify a missassembly from an alignment file.

    # Args
    * `aln`
        * Alignment file as BAM or CRAM file. Requires fasta and `cs` tag if CRAM.
    * `bed`
        * BED3 file with regions to evaluate.
    * `ignore_bed`
        * BED3 file with regions to ignore.
    * `threads`
        * Number of threads to spawn.
    * `cfg`
        * Configfile. See [`nucflag::config::Config`]
    * `preset`
        * Configuration for specific LR sequencing reads.
        * Modifies `cfg` where preset specific options take priority.
        * See [`nucflag::preset::Preset`].

    # Returns
    * [`PyNucFlagResult`]
    """

def run_nucflag_itv(
    aln: str,
    itv: tuple[int, int, str],
    fasta: str | None = None,
    ignore_bed: str | None = None,
    threads: int = 1,
    cfg: str | None = None,
    preset: str | None = None,
) -> PyNucFlagResult:
    """Classify a missassembly for one interval. Identical to `run_nucflag` but only for one interval."""
