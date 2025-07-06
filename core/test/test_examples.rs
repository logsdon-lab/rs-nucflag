use coitrees::Interval;
use polars::{
    io::SerReader,
    prelude::{CsvParseOptions, CsvReadOptions},
};
use rs_nucflag::{
    io::{read_cfg, write_tsv},
    nucflag,
};

fn check_output(
    aln: &str,
    bed: &str,
    fasta: Option<&str>,
    config: Option<&str>,
    expected: Option<&str>,
    save_res: Option<&str>,
) {
    let itvs = std::fs::read_to_string(bed).unwrap();
    let itv = itvs
        .lines()
        .next()
        .map(|itv| {
            let [chrom, st, end] = itv.split('\t').collect::<Vec<&str>>()[..] else {
                panic!("Invalid bed.")
            };
            Interval::new(
                st.parse::<i32>().unwrap(),
                end.parse::<i32>().unwrap(),
                chrom.to_owned(),
            )
        })
        .unwrap();
    let mut res = nucflag(
        aln,
        fasta,
        &itv,
        None,
        config
            .map(|cfg| read_cfg(Some(cfg), None).unwrap())
            .unwrap_or_default(),
    )
    .unwrap();

    if save_res.is_some() {
        write_tsv(&mut res.regions, save_res).unwrap();
    }

    if let Some(expected) = expected {
        let df_expected = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
            .try_into_reader_with_file_path(Some(expected.into()))
            .unwrap()
            .finish()
            .unwrap();

        assert_eq!(
            res.regions, df_expected,
            "Called regions for ({aln}) not equal."
        );
    }
}

#[test]
fn test_dupes() {
    let indir = "test/dupes/input";
    let expdir = "test/dupes/expected";
    for case in ["aln_1", "aln_2", "aln_3"] {
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, None, None, Some(&expected));
        check_output(&aln, &bed, None, None, Some(&expected), None)
    }
}

#[test]
fn test_ending_scaffold() {
    let indir = "test/ending_scaffold/input";
    let expdir = "test/ending_scaffold/expected";
    {
        let case = "aln_1";
        let aln = format!("{indir}/{case}.bam");
        let fa = format!("{indir}/{case}.fa");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, Some(&fa), None, None, Some(&expected));
        check_output(&aln, &bed, Some(&fa), None, Some(&expected), None)
    }
}

#[test]
fn test_ignore_low_cov_boundaries() {
    let indir = "test/ignore_boundaries/input";
    let expdir = "test/ignore_boundaries/expected";
    for (case, suffix) in [("aln_1", "_ignored"), ("aln_1", "")] {
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let cfg = format!("{indir}/{case}{suffix}.toml");
        let expected = format!("{expdir}{suffix}/{case}.bed");
        // check_output(&aln, &bed, None, Some(&cfg), None, Some(&expected));
        check_output(&aln, &bed, None, Some(&cfg), Some(&expected), None)
    }
}

#[test]
fn test_het() {
    let indir = "test/het/input";
    let expdir = "test/het/expected";
    {
        let case = "aln_1";
        let aln = format!("{indir}/{case}.bam");
        let cfg = format!("{indir}/{case}.toml");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, Some(&cfg), None, Some(&expected));
        check_output(&aln, &bed, None, Some(&cfg), Some(&expected), None)
    }
}

#[test]
fn test_hsat() {
    let indir = "test/hsat/input";
    let expdir = "test/hsat/expected";
    {
        let case = "aln_1";
        let aln = format!("{indir}/{case}.bam");
        let fa = format!("{indir}/{case}.fa");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, Some(&fa), None, None, Some(&expected));
        check_output(&aln, &bed, Some(&fa), None, Some(&expected), None)
    }
}

#[test]
fn test_minor_collapse() {
    let indir = "test/minor_collapse/input";
    let expdir = "test/minor_collapse/expected";
    for case in ["aln_1", "aln_2", "aln_3", "aln_4"] {
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, None, None, Some(&expected));
        check_output(&aln, &bed, None, None, Some(&expected), None)
    }
}

#[test]
fn test_misjoin() {
    let indir = "test/misjoin/input";
    let expdir = "test/misjoin/expected";
    {
        let case = "aln_1";
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, None, None, Some(&expected));
        check_output(&aln, &bed, None, None, Some(&expected), None)
    }
}

#[test]
fn test_standard() {
    let indir = "test/standard/input";
    let expdir = "test/standard/expected";
    {
        let case = "aln_1";
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, None, None, Some(&expected));
        check_output(&aln, &bed, None, None, Some(&expected), None)
    }
}

#[test]
fn test_ignore_false_collapse() {
    let indir = "test/ignore_false_collapse/input";
    let expdir = "test/ignore_false_collapse/expected";
    for case in ["aln_1", "aln_2"] {
        let aln = format!("{indir}/{case}.cram");
        let fa = format!("{indir}/{case}.fa");
        let cfg = format!("{indir}/{case}.toml");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, Some(&fa), Some(&cfg), None, Some(&expected));
        check_output(&aln, &bed, Some(&fa), Some(&cfg), Some(&expected), None)
    }
}
