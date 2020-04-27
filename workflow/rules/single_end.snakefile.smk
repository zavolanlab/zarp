rule remove_adapters_cutadapt:
    '''
        Remove adapters
    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "start",
            "{sample}.fq1.fastq.gz")

    output:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_adapters_mate1.fastq.gz")

    params:
        adapters_3 = lambda wildcards:
            samples_table.loc[wildcards.sample, 'fq1_3p'],
        adapters_5 = lambda wildcards:
            samples_table.loc[wildcards.sample, 'fq1_5p']

    singularity:
        "docker://zavolab/cutadapt:1.16-slim"

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_adapters_cutadapt.se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_adapters_cutadapt.se.stdout.log")
    shell:
        "(cutadapt \
        -e 0.1 \
        -j {threads} \
        -m 10 \
        -n 2 \
        -a {params.adapters_3} \
        -g {params.adapters_5} \
        -o {output.reads} \
        {input.reads}) \
        1> {log.stdout} 2> {log.stderr}"


rule remove_polya_cutadapt:
    '''
        Remove ployA  tails
    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_adapters_mate1.fastq.gz")

    output:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_polya_mate1.fastq.gz")

    params:
        polya_3 = lambda wildcards:
            samples_table.loc[wildcards.sample, "fq1_polya_3p"],
        polya_5 = lambda wildcards:
            samples_table.loc[wildcards.sample, "fq1_polya_5p"]

    singularity:
        "docker://zavolab/cutadapt:1.16-slim"

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_polya_cutadapt.se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_polya_cutadapt.se.stdout.log")

    shell:
        "(cutadapt \
        -j {threads} \
        -n 1 \
        -e 0.1 \
        -O 1 \
        -m 10  \
        -a {params.polya_3} \
        -g {params.polya_5} \
        -o {output.reads} \
        {input.reads};) \
        1> {log.stdout} 2> {log.stderr}"


rule map_genome_star:
    '''
        Map to genome using STAR
    '''
    input:
        index = lambda wildcards:
            os.path.join(
                config["star_indexes"],
                str(samples_table.loc[wildcards.sample, "organism"]),
                str(samples_table.loc[wildcards.sample, "index_size"]),
                "STAR_index",
                "chrNameLength.txt"),
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_polya_mate1.fastq.gz")

    output:
        bam = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.se.Aligned.sortedByCoord.out.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.se.Log.final.out")

    params:
        sample_id = "{sample}",
        index = lambda wildcards:
            os.path.join(
                config["star_indexes"],
                str(samples_table.loc[wildcards.sample, "organism"]),
                str(samples_table.loc[wildcards.sample, "index_size"]),
                "STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.se."),
        multimappers = lambda wildcards:
                samples_table.loc[wildcards.sample, "multimappers"],
        soft_clip = lambda wildcards:
                samples_table.loc[wildcards.sample, "soft_clip"],
        pass_mode = lambda wildcards:
                samples_table.loc[wildcards.sample, "pass_mode"],

    singularity:
        "docker://zavolab/star:2.7.3a-slim"

    threads: 12

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "map_genome_star.se.stderr.log")

    shell:
        "(STAR \
        --runMode alignReads \
        -- twopassMode {params.pass_mode} \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.reads} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 1 \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_SortedByCoordinate \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        --outFilterType BySJout \
        --outReadsUnmapped None \
        --outSAMattrRGline ID:rnaseq_pipeline SM:{params.sample_id} \
        --alignEndsType {params.soft_clip} > {output.bam};) \
        2> {log.stderr}"


rule quantification_salmon:
    '''
        Quantification at transcript and gene level using Salmon
    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_polya_mate1.fastq.gz"),
        index = lambda wildcards:
            os.path.join(
                config["salmon_indexes"],
                str(samples_table.loc[wildcards.sample, "organism"]),
                str(samples_table.loc[wildcards.sample, "kmer"]),
                "salmon.idx"),
        gtf = lambda wildcards:
            samples_table.loc[wildcards.sample, "gtf"]

    output:
        gn_estimates = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se",
            "quant.genes.sf"),
        tr_estimates = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se",
            "quant.sf")

    params:
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se"),
        libType = lambda wildcards:
                samples_table.loc[wildcards.sample, "libtype"]

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "quantification_salmon.se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "quantification_salmon.se.stdout.log")

    threads: 12

    singularity:
        "docker://zavolab/salmon:1.1.0-slim"

    shell:
        "(salmon quant \
        --libType {params.libType} \
        --seqBias \
        --validateMappings \
        --threads {threads} \
        --writeUnmappedNames \
        --index {input.index} \
        --geneMap {input.gtf} \
        --unmatedReads {input.reads} \
        -o {params.output_dir};) \
        1> {log.stdout} 2> {log.stderr}"


rule genome_quantification_kallisto:
    '''
        Quantification at transcript and gene level using Kallisto
    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_polya_mate1.fastq.gz"),
        index = lambda wildcards:
            os.path.join(
                config["kallisto_indexes"],
                samples_table.loc[wildcards.sample, "organism"],
                "kallisto.idx")

    output:
        pseudoalignment = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto",
            "{sample}.se.kallisto.pseudo.sam")

    params:
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto"),
        fraglen = lambda wildcards:
            samples_table.loc[wildcards.sample, 'mean'],
        fragsd = lambda wildcards:
            samples_table.loc[wildcards.sample, 'sd'],
        directionality = lambda wildcards:
            samples_table.loc[wildcards.sample, 'kallisto_directionality']

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "genome_quantification_kallisto.se.stderr.log")

    singularity:
        "docker://zavolab/kallisto:0.46.1-slim"

    shell:
        "(kallisto quant \
        -i {input.index} \
        -o {params.output_dir} \
        --single \
        -l {params.fraglen} \
        -s {params.fragsd} \
        --pseudobam \
        {params.directionality}-stranded \
        {input.reads} > {output.pseudoalignment};) \
        2> {log.stderr}"

