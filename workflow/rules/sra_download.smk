"""Independent Snakefile for downloading samples from SRA."""

import pandas as pd

samples = pd.read_csv(config["samples"], header=0, index_col=0, sep = "\t")
DOWNLOAD_DIR = config["outdir"]
# Write fastq.gz location into new sample table.
SAMPLES_OUT = config["samples_out"]
samples_mod = samples.copy()

localrules: prefetch, add_fq_file_path, all

rule all:
    "Target rule."
    input:
        SAMPLES_OUT


rule prefetch:
    "Prefetch SRA entry. Requires internet access."
    output:
        os.path.join(DOWNLOAD_DIR, "{sample}", "{sample}.sra")
    params:
        outdir = DOWNLOAD_DIR
    conda:
        os.path.join(workflow.basedir, "..", "envs", "sra-tools.yaml")
    singularity:
        "docker://ncbi/sra-tools"
    shell:
        """
        prefetch {wildcards.sample} --output-directory {params.outdir}
        """

rule fasterq_dump:
    "Dump SRA entry as fastq file(s)."
    input:
        os.path.join(DOWNLOAD_DIR, "{sample}", "{sample}.sra")
    output:
        os.path.join(DOWNLOAD_DIR, "{sample}", "{sample}.dumped")
    params:
        outdir = lambda wildcards: os.path.join(DOWNLOAD_DIR, wildcards.sample),
        cluster_log_path = config["cluster_log_dir"]
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    threads: 4
    conda:
        os.path.join(workflow.basedir, "..", "envs", "sra-tools.yaml")
    singularity:
        "docker://ncbi/sra-tools"
    shell:
        """
        fasterq-dump {params.outdir} --outdir {params.outdir} \
            --mem {resources.mem_mb}MB --threads {threads} \
            --temp {resources.tmpdir}; \
        touch {output}
        """


def get_fastq_files(wildcards):
    """Obtain paths to fastq files for given sample.

    Args:
        wildcards: Snakemake wildcards, contains name "sample".
    
    Returns:
        list (str): paths to .fastq files.
    """
    files = os.listdir(os.path.join(DOWNLOAD_DIR, wildcards.sample))
    to_zip = []
    for f in files:
        if f.endswith(".fastq"):
            to_zip.append(os.path.join(DOWNLOAD_DIR, wildcards.sample, f))
    return(to_zip)


rule compress_fastq:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        tmpf = os.path.join(DOWNLOAD_DIR, "{sample}", "{sample}.dumped")
    output:
         os.path.join(DOWNLOAD_DIR, "{sample}", "{sample}.processed")
    params:
        files = get_fastq_files,
        cluster_log_path = config["cluster_log_dir"]
    threads: 6
    conda:
        os.path.join(workflow.basedir, "..", "envs", "pigz.yaml")
    singularity:
        "docker://bytesco/pigz"
    shell:
        """
        pigz --best --processes {threads} {params.files}; \
        touch {output}
        """

rule add_fq_file_path:
    "Add fastq paths to sample table."
    input:
        expand(os.path.join(DOWNLOAD_DIR, 
            "{sample}", "{sample}.processed"), 
            sample = samples[samples.index.str.contains("SRR")].index.tolist())
    output:
        SAMPLES_OUT
    run:
        for sample in input:
            files = os.listdir(os.path.dirname(sample))
            sample_name = os.path.basename(sample).split(".")[0]
            gzs = []
            for f in files:
                if f.endswith(".fastq.gz"):
                    gzs.append(os.path.join(DOWNLOAD_DIR, sample_name, f))
            if len(gzs) == 1:
                # single-end sample
                samples_mod.loc[sample_name, "fq1"] = gzs[0]
            if len(gzs) == 2:
                gzs.sort()
                samples_mod.loc[sample_name, "fq1"] = gzs[0]
                samples_mod.loc[sample_name, "fq2"] = gzs[1]
        samples_mod.to_csv(SAMPLES_OUT, index=True, sep = "\t")
