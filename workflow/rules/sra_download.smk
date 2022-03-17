"""Independent Snakefile for downloading samples from SRA."""

localrules: prefetch, all

rule prefetch:
    "Requires internet access."
    output:
        os.path.join("sra_downloads", "{sample}", "{sample}.sra")
    params:
        outdir = "sra_downloads"
    conda:
        "../envs/sra-tools.yaml"
    singularity:
        "docker://ncbi/sra-tools"
    shell:
        """
        prefetch {wildcards.sample} --output-directory {params.outdir}
        """

rule fasterq_dump:
    input:
        os.path.join("sra_downloads", "{sample}", "{sample}.sra")
    output:
        #directory(os.path.join("sra_downloads", "{sample}"))
        os.path.join("sra_downloads", "{sample}", "{sample}.dumped")
    params:
        outdir = lambda wildcards: os.path.join("sra_downloads", wildcards.sample)
    resources:
        mem_mb = 2000
    threads: 6
    conda:
        "../envs/sra-tools.yaml"
    singularity:
        "docker://ncbi/sra-tools"
    shell:
        """
        fasterq-dump {params.outdir} --outdir {params.outdir} \
            --mem {resources.mem_mb}MB --threads {threads} \
            --temp {resources.tmpdir}; \
        touch {output}
        """

# def get_fastq_files(wildcards):
#     files = os.listdir(checkpoints.fasterq_dump.get(**wildcards).output[0])
#     to_zip = []
#     for f in files:
#         if f.endswith(".fastq"):
#             to_zip.append(os.path.join("sra_downloads", wildcards.sample, f))
#     if len(to_zip) == 1:
#         return expand(os.path.join("sra_downloads", "{sample}", "{sample}.fastq"),
#             sample = wildcards.sample)
#     elif len(to_zip) == 2: 
#         return expand(os.path.join("sra_downloads", "{sample}", "{sample}_{mate}.fastq"),
#             sample = wildcards.sample,
#             mate = ["1", "2"])

def get_fastq_files(wildcards):
    files = os.listdir(os.path.join("sra_downloads", wildcards.sample))
    to_zip = []
    for f in files:
        if f.endswith(".fastq"):
            to_zip.append(os.path.join("sra_downloads", wildcards.sample, f))
    return(to_zip)


rule compress_fastq:
    "Compress with pigz at best (9) compression level."
    input:
        # files = glob_wildcards(os.path.join("sra_downloads", "{sample}", "{id}.fastq"))
        tmpf = os.path.join("sra_downloads", "{sample}", "{sample}.dumped")
        # os.path.join("sra_downloads", "{sample}", "{sample}{id}.fastq")
    output:
         os.path.join("sra_downloads", "{sample}", "{sample}.processed")
    params:
        files = get_fastq_files
    threads: 6
    conda:
        "../envs/pigz.yaml"
    singularity:
        "docker://bytesco/pigz"
    shell:
        """
        pigz --best --processes {threads} {params.files}; \
        touch {output}
        """

rule all:
    input:
        expand(os.path.join("sra_downloads", "{sample}", "{sample}.processed"), 
            sample = ["SRR179707", "SRR2969253"])