rule fasterq_dump_pe:
    "Dump SRA entry as fastq file(s)."
    input:
        infile=os.path.join(config["outdir"], "prefetch", "{sample}"),
    output:
        outfile1=temp(
            os.path.join(
                config["outdir"], "fasterq_dump", "{sample}", "{sample}_1.fastq"
            )
        ),
        outfile2=temp(
            os.path.join(
                config["outdir"], "fasterq_dump", "{sample}", "{sample}_2.fastq"
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
    threads: 4
    conda:
        os.path.join(workflow.basedir, "..", "envs", "sra-tools.yaml")
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "fasterq_dump.pe.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "fasterq_dump.pe.stdout.log"
        ),
    shell:
        """
        (mkdir -p {params.sample_dir};\
        cd {params.results}; \
        fasterq-dump {params.sample} \
        --mem {resources.mem_mb}MB \
        --outdir {params.outdir} \
        --threads {threads} \
        --temp tmpdir;) 1> {log.stdout} 2> {log.stderr}
        """


rule compress_fastq_pe:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        file1=os.path.join(
            config["outdir"], "fasterq_dump", "{sample}", "{sample}_1.fastq"
        ),
        file2=os.path.join(
            config["outdir"], "fasterq_dump", "{sample}", "{sample}_2.fastq"
        ),
    output:
        file1=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}_1.fastq.gz"
        ),
        file2=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}_2.fastq.gz"
        ),
    params:
        cluster_log_path=config["cluster_log_dir"],
        outdir=os.path.join(config["outdir"], "compress", "{sample}"),
    threads: 6
    conda:
        os.path.join(workflow.basedir, "..", "envs", "pigz.yaml")
    container:
        "docker://quay.io/biocontainers/pigz:2.8"
    log:
        stderr=os.path.join(
            config["log_dir"], "samples", "{sample}", "compress__pe_fastq.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "samples", "{sample}", "compress_pe_fastq.stdout.log"
        ),
    shell:
        """ (mkdir -p {params.outdir}; \
            pigz --best --processes {threads} {input.file1} --stdout > {output.file1}; \
            pigz --best --processes {threads} {input.file2} --stdout > {output.file2};) \
            1> {log.stdout} 2> {log.stderr};
        """


rule process_fastq_pe:
    "Keep the fastq.gz file paths in a table"
    input:
        file1=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}_1.fastq.gz"
        ),
        file2=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}_2.fastq.gz"
        ),
    output:
        outfile=os.path.join(
            config["outdir"], "compress", "{sample}", "{sample}.pe.tsv"
        ),
    params:
        cluster_log_path=config["cluster_log_dir"],
        filename="{sample}",
    threads: 1
    log:
        stderr=os.path.join(
            config["log_dir"], "{sample}", "process_pe_fastq.stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"], "{sample}", "process_pe_fastq.stdout.log"
        ),
    run:
        samples_mod = pd.DataFrame()
        samples_mod.index.name = "sample"
        samples_mod["fq1"] = ""
        samples_mod["fq2"] = ""
        samples_mod.loc[params.filename, "fq1"] = input.file1
        samples_mod.loc[params.filename, "fq2"] = input.file2
        samples_mod.to_csv(output.outfile, index=True, sep="\t")
