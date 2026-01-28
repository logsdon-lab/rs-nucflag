cargo run --release --manifest-path examples/Cargo.toml -- \
    core/test/collapse/input/aln_1.bam \
    core/test/collapse/input/aln_1.fa.gz \
    core/test/collapse/input/aln_1.bed \
    core/nucflag.toml \
    1 \
    none
