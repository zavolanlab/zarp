"""Independent Snakefile for downloading samples from SRA."""

import pandas as pd


include: os.path.join("sra_paired_end.smk")
include: os.path.join("sra_single_end.smk")


samples = pd.read_csv(config["samples"], header=0, index_col=0, sep="\t")


localrules:
    all,
    prefetch,
    add_fq_file_path,
    get_layout,


rule all:
    "Target rule."
    input:
        final=config["samples_out"],


checkpoint get_layout:
    "Get the library type of each sample (paired or single-end)."
    output:
        outdir=directory(os.path.join(config["outdir"], "get_layout", "{sample}")),
    params:
        cluster_log_path=config["cluster_log_dir"],
    conda:
        os.path.join(workflow.basedir, "..", "envs", "entrez-direct.yaml")
    singularity:
        "docker://quay.io/biocontainers/entrez-direct:16.2--he881be0_1"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "get_layout.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "get_layout.stdout.log"
        ),
    shell:
        """
        (mkdir -p {output.outdir}; \
        layout=$(efetch -db sra \
        -id {wildcards.sample} \
        -format runinfo | \
        csv2xml -set Set -rec Rec -header | \
        xtract -pattern Rec -if Run -equals {wildcards.sample} -element LibraryLayout); \
        touch {output.outdir}/$layout.info ; \
        ) 1> {log.stdout} 2> {log.stderr}
        """


rule prefetch:
    "Prefetch SRA entry. Requires internet access."
    output:
        outdir=directory(os.path.join(config["outdir"], "prefetch", "{sample}")),
    params:
        cluster_log_path=config["cluster_log_dir"],
        outdir=os.path.join(config["outdir"], "prefetch"),
    conda:
        os.path.join(workflow.basedir, "..", "envs", "sra-tools.yaml")
    singularity:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "prefetch.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "prefetch.stdout.log"
        ),
    shell:
        """
        (mkdir -p {params.outdir}; \
        cd {params.outdir}; \
        prefetch {wildcards.sample} \
        ) 1> {log.stdout} 2> {log.stderr}
        """


def get_layouts(wildcards):
    ivals = []
    for i in samples[
        samples.index.str.contains("^.RR", regex=True, case=True)
    ].index.tolist():
        checkpoint_output = checkpoints.get_layout.get(
            sample=i, **wildcards
        ).output.outdir
        ivals.extend(
            glob_wildcards(os.path.join(checkpoint_output, "{layout}.info")).layout
        )
    ivals2 = []
    for ival in ivals:
        if ival == "PAIRED":
            ivals2.append("pe")
        elif ival == "SINGLE":
            ivals2.append("se")

    layouts = expand(
        os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}.{seqmode}.tsv"
        ),
        zip,
        sample=[
            i
            for i in samples[
                samples.index.str.contains("^.RR", regex=True, case=True)
            ].index.tolist()
        ],
        seqmode=ivals2,
    )
    return layouts


rule add_fq_file_path:
    "Add fastq paths to sample table."
    input:
        files=get_layouts,
    output:
        outfile=config["samples_out"],
    params:
        cluster_log_path=config["cluster_log_dir"],
    run:
        samples_mod = []
        for eachtable in input.files:
            table = pd.read_csv(eachtable, sep="\t", header=0, index_col=0)
            samples_mod.append(table)
        samples_mod = pd.concat(samples_mod)
        samples_mod.to_csv(output.outfile, index=True, sep="\t")
