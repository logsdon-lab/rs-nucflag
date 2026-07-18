# `rs-nucflag`
Rust implementation of `NucFlag`

## Structure
This repo is structured as follows:
```
core/
├── Cargo.lock
├── Cargo.toml
├── LICENSE
├── nucflag.toml                # Example configuration file
├── README.md
├── src
│   ├── binning.rs              # Partitioning of pileup
│   ├── classify.rs             # Classification of calls
│   ├── config.rs               # Configuration structs
│   ├── intervals.rs            # Interval helpers
│   ├── io.rs                   # Input/output helpers
│   ├── lib.rs                  # Module imports
│   ├── misassembly.rs          # Misassembly types
│   ├── nucflag.rs              # Entrypoint
│   ├── peak.rs                 # Peak calling
│   ├── pileup.rs               # Pileup
│   ├── postprocess.rs          # Call postprocessing and merging
│   ├── preset.rs               # Configuration presets
│   └── repeats.rs              # Repeat detection
└── test
    ├── collapse                # Example of a collapse
    ├── dupes                   # Example of a false duplication
    ├── ending_scaffold         # Example of omission of an ending scaffold
    ├── het                     # Example of a heterozygous site
    ├── hsat                    # Example of HSAT-1A coverage bias
    ├── ignore_boundaries       # Example of ignoring boundary errors
    ├── ignore_false_collapse   # Example of ignoring a false collapse due to region boundary
    ├── inversion               # Example of a false inversion
    ├── minor_collapse          # Example of a small collapse
    ├── misjoin                 # Example of a misjoin
    ├── no_indel_merging        # Example to show that indels should not be merged
    ├── pileup                  # Test files for pileup unit test
    ├── regenerate_plots.json   # JSON manifest of files for workflow to regenerate plots
    ├── regenerate_plots.smk    # Workflow to regenerate plots
    ├── standard                # Example used in README
    └── test_examples.rs        # All tests
```

Each `test` subdirectory contains:
* `input`
    * BAM/CRAM/BED/TOML files for tests
* `expected`
    * Output plots and BED files.

## Commands
Remake NucFreq plots. Requires `snakemake` and `nucflag>=1.0`
```bash
which snakemake nucflag
make test_rs_remake_images
```

Run all tests.
```bash
make test_rs
```

Build library.
```bash
make build_rs
```
