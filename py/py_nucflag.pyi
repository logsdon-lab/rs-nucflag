import polars as pl

# flake8: noqa: PYI021

class PyNucFlagResult:
    """
    `NucFlag` results.
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
    Pileup of region with columns:
    * `chrom`
        * Chromosome name
    * `pos`
        * Position
    * `cov`
        * Coverage
    * `status`
        * Status
    * `mismatch`
        * Number of mismatches
    * `mapq`
        * MAPQ
    * `insertion`
        * Number of insertions
    * `deletion`
        * Number of deletions
    * `softclip`
        * Number of softclipped bases
    * `bin`
        * Region number.
        * 0 if non-repetitive and some number if repetitive
    * `bin_ident`
        * Self-sequence identity of region.
        * 0.0 if non-repetitive
    """
    regions: pl.DataFrame
    """
    Regions and their status in BED9 format:
    * `#chrom`
        * Chromosome name
    * `chromStart`
        * Chromosome start
    * `chromEnd`
        * Chromosome end
    * `name`
        * Name of status.
    * `score`
        * Coverage of region.
    * `strand`
        * None or `.`
    * `thickStart`
        * Same as `chromEnd`
    * `thickEnd`
        * Same as `chromStart`
    * `itemRgb`
        * Color of status
    * `zscore`
        * [MAD adjusted z-score.](https://www.ibm.com/docs/en/cognos-analytics/12.0.x?topic=terms-modified-z-score)
    * `af`
        * Allele frequency.
        * `signal / coverage` where deletions are not counted toward `coverage`
    """

def get_regions(
    aln: str, bed: str | None = None, window: int = 10000000
) -> list[tuple[int, int, str]]:
    """
    Get interval regions from an alignment file or bed file.

    # Args
    * `aln`
        * Alignment as BAM or CRAM file.
    * `bed`
        * BED file with coordinates.
        * If not provided, splits all regions listed in `aln` header (`@SQ`) into non-overlapping windows.
        * Invalid intervals are ignored.
    * `window`
        * Window size in base pairs of non-overlapping windows if `bed` not provided.

    # Returns
    * List of intervals as tuples with the format, start, end, and chromosome name.

    # Example
    ```python
    from py_nucflag import get_regions

    # Get all possible aligned regions
    regions = get_regions("aln.bam", window = 10_000_000)
    # Get regions in bed file and verify that exist from alignment header.
    regions_bed = get_regions("aln.bam", bed="regions.bed")
    ```
    """

def get_config_from_preset(preset: str | None = None, cfg: str | None = None) -> str:
    """
    Return `NucFlag` config as TOML string from preset.

    # Args
    * `preset`
        * `NucFlag` preset.
        * Either `ont_r9`, `ont_r10`, or `hifi`
    * `config`
        * Path to configfile.
        * If provided alongside `preset`, `preset` fields given priority.

    # Returns
    * Configuration as a TOML string.

    # Example
    ```python
    from py_nucflag import get_config_from_preset

    # Get preset parameters
    config = get_config_from_preset(preset="hifi")

    # Read in existing config and merge with preset parameters.
    config_w_preset = get_config_from_preset(preset="ont_r10", cfg="config.toml")
    ```
    """

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

    # Example
    ```python
    from py_nucflag import run_nucflag

    # Run on all alignments.
    # Providing a bed file is highly recommended as pileup information is stored.
    results = run_nucflag(
        "sample.bam",
        bed="regions.bed",
        fasta="sample.fa.gz"
    )
    # Iterate through each region's misassembly calls.
    for result in results:
        # The misassembly calls as BED9
        print(result.regions)
        # The raw pileup
        print(result.pileup)
    ```
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
    """
    Classify a missassembly for one interval. Identical to `run_nucflag` but only for one interval.

    # Example
    ```python
    from py_nucflag import run_nucflag_itv

    # Run on alignments to a single interval.
    result = run_nucflag_itv(
        "sample.bam",
        itv=(0, 1_000_000), "chr1"),
        fasta="sample.fa.gz"
    )
    print(result.regions)
    ```
    """
