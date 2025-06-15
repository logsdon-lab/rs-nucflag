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
    for case in ["aln_1"] {
        let aln = format!("{indir}/{case}.bam");
        let fa = format!("{indir}/{case}.fa");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, Some(&fa), None, None, Some(&expected));
        check_output(&aln, &bed, Some(&fa), None, Some(&expected), None)
    }
}

#[test]
fn test_het() {
    let indir = "test/het/input";
    let expdir = "test/het/expected";
    for case in ["aln_1"] {
        let aln = format!("{indir}/{case}.bam");
        let cfg = format!("{indir}/{case}.toml");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, Some(&cfg), None, Some(&expected));
        check_output(&aln, &bed, None, Some(&cfg), Some(&expected), None)
    }
}

// #[test]
// fn test_hsat() {
//     todo!("Figure out what is nondetermenistic in output for good coverage");
//     /*
//     diff /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/hsat/expected/aln_1.bed /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/hsat/expected/aln_1_.bed
//     25c25
//     < NA18534_chr4_haplotype2-0000075       52388511        53432829        good    19.0    +       52388511        53432829        0,0,0
//     ---
//     > NA18534_chr4_haplotype2-0000075       52388511        53432829        good    21.0    +       52388511        53432829        0,0,0
//     */
//     let indir = "test/hsat/input";
//     let expdir = "test/hsat/expected";
//     for case in ["aln_1"] {
//         let aln = format!("{indir}/{case}.bam");
//         let fa = format!("{indir}/{case}.fa");
//         let bed = format!("{indir}/{case}.bed");
//         let expected = format!("{expdir}/{case}.bed");
//         // check_output(&aln, &bed, Some(&fa), None, None, Some(&expected));
//         check_output(&aln, &bed, Some(&fa), None, Some(&expected), Some(&format!("{expdir}/{case}_.bed")))
//     }
// }

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
    for case in ["aln_1"] {
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
    for case in ["aln_1"] {
        let aln = format!("{indir}/{case}.bam");
        let bed = format!("{indir}/{case}.bed");
        let expected = format!("{expdir}/{case}.bed");
        // check_output(&aln, &bed, None, None, None, Some(&expected));
        check_output(&aln, &bed, None, None, Some(&expected), None)
    }
}
