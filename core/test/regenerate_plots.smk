import os
import json
from os.path import join


with open("regenerate_plots.json", "rt") as fh:
    cases = json.load(fh)


rule run_nucflag:
    input:
        bam=lambda wc: cases[wc.case]["bam"],
        bed=lambda wc: cases[wc.case]["bed"],
        fasta=lambda wc: (
            cases[wc.case]["fasta"]
            if cases[wc.case].get("fasta", [])
            else []
        ),
        config=lambda wc: (
            cases[wc.case]["config"]
            if cases[wc.case].get("config", [])
            else []
        ),
    output:
        touch(join("results", "{case}.done"))
    resources:
        mem="50GB"
    params:
        fasta=lambda wc, input: f"-f {input.fasta}" if input.fasta else "",
        config=lambda wc, input: f"-c {input.config}" if input.config else "",
        output_dir=lambda wc: os.path.dirname(cases[wc.case]["output_png"]),
        output_png=lambda wc: cases[wc.case]["output_png"],
        final_png=lambda wc: cases[wc.case]["final_png"],
        additional_flags=lambda wc: " ".join(cases[wc.case]["additional_flags"]),
    threads: 4
    shell:
        """
        nucflag \
        -i {input.bam} \
        -b {input.bed} \
        -d {params.output_dir} \
        -t {threads} \
        -p {threads} {params.config} {params.fasta} {params.additional_flags}
        mv {params.output_png} {params.final_png}
        """

rule all:
    input:
        expand(
            rules.run_nucflag.output,
            case=cases.keys(),
        ),
    default_target: True
