rule fasterq_dump_se:
    "Dump SRA entry as fastq file(s)."
    input:
        infile=os.path.join(config["outdir"], "prefetch", "{sample}"),
    output:
        outfile1=temp(
            os.path.join(
                config["outdir"], "fasterq_dump", "{sample}", "{sample}.fastq"
            )
        ),
    params:
        cluster_log_path=config["cluster_log_dir"],
        results=os.path.join(config["outdir"], "fasterq_dump"),
        outdir="{sample}",
        sample_dir=os.path.join(config["outdir"], "fasterq_dump", "{sample}"),
        sample=os.path.abspath(os.path.join(config["outdir"], "prefetch", "{sample}")),
    resources:
        mem_mb=lambda wildcards, attempt: 3048 * attempt,
        tmpdir=os.path.join("tmpdir"),
    threads: 4
    conda:
        os.path.join(workflow.basedir, "..", "envs", "sra-tools.yaml")
    singularity:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "fasterq_dump.se.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "fasterq_dump.se.stdout.log"
        ),
    shell:
        """
        (mkdir -p {params.sample_dir};\
        cd {params.results}; \
        fasterq-dump {params.sample} \
        --mem {resources.mem_mb}MB \
        --outdir {params.outdir} \
        --threads {threads} \
        --temp {resources.tmpdir};) 1> {log.stdout} 2> {log.stderr}
        """


rule compress_fastq_se:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        file1=os.path.join(
            config["outdir"], "fasterq_dump", "{sample}", "{sample}.fastq"
        ),
    output:
        file1=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}.fastq.gz"
        ),
    params:
        cluster_log_path=config["cluster_log_dir"],
        outdir=os.path.join(config["outdir"], "compress", "{sample}"),
    threads: 6
    conda:
        os.path.join(workflow.basedir, "..", "envs", "pigz.yaml")
    singularity:
        "docker://quay.io/biocontainers/pigz:2.8"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "compress_se_fastq.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "compress_se_fastq.stdout.log"
        ),
    shell:
        """(mkdir -p {params.outdir}; \
        pigz --best --processes {threads} {input.file1} --stdout > {output.file1};) \
        1> {log.stdout} 2> {log.stderr};
        """


rule process_fastq_se:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        file=os.path.join(config["outdir"], "compress", "{sample}", "{sample}.fastq.gz"),
    output:
        outfile=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}.se.tsv"
        ),
    params:
        cluster_log_path=config["cluster_log_dir"],
        filename="{sample}",
    threads: 6
    log:
        stderr=os.path.join(config["log_dir"], "{sample}__process_se_fastq.stderr.log"),
        stdout=os.path.join(config["log_dir"], "{sample}__process_se_fastq.stdout.log"),
    run:
        samples_mod = pd.DataFrame()
        samples_mod.index.name = "sample"
        samples_mod["fq1"] = ""
        samples_mod["fq2"] = ""
        samples_mod.loc[params.filename, "fq1"] = input.file
        samples_mod.to_csv(output.outfile, index=True, sep="\t")
