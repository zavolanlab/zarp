rule pe_remove_adapters_cutadapt:
    '''
        Remove adapters
    '''
    input:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "start",
            "{sample}.fq1.fastq.gz"),

        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "start",
            "{sample}.fq2.fastq.gz"),

    output:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate2.fastq.gz")

    params:
        adapter_3_mate1 = lambda wildcards:
            get_sample('fq1_3p', search_id='index', search_value=wildcards.sample),
        adapter_5_mate1 = lambda wildcards:
            get_sample('fq1_5p', search_id='index', search_value=wildcards.sample),
        adapter_3_mate2 = lambda wildcards:
            get_sample('fq2_3p', search_id='index', search_value=wildcards.sample),
        adapter_5_mate2 = lambda wildcards:
            get_sample('fq2_5p', search_id='index', search_value=wildcards.sample)

    singularity:
        "docker://zavolab/cutadapt:1.16-slim"

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_adapters_cutadapt.pe.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_adapters_cutadapt.pe.stdout.log")

    shell:
        "(cutadapt \
        -e 0.1 \
        -j {threads} \
        --pair-filter=any \
        -m 10 \
        -n 2 \
        -a {params.adapter_3_mate1} \
        -g {params.adapter_5_mate1} \
        -A {params.adapter_3_mate2} \
        -G {params.adapter_5_mate2} \
        -o {output.reads1} \
        -p {output.reads2} \
        {input.reads1} \
        {input.reads2};) \
        1> {log.stdout} 2>{log.stderr}"


rule pe_remove_polya_cutadapt:
    '''
        Remove polyA tails
    '''
    input:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate2.fastq.gz")

    output:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate2.fastq.gz")

    params:
        polya_3_mate1 = lambda wildcards:
            get_sample(
                'fq1_polya_3p',
                search_id='index',
                search_value=wildcards.sample),
        polya_5_mate1 = lambda wildcards:
            get_sample(
                'fq1_polya_5p',
                search_id='index',
                search_value=wildcards.sample),
        polya_3_mate2 = lambda wildcards:
            get_sample(
                'fq2_polya_3p',
                search_id='index',
                search_value=wildcards.sample),
        polya_5_mate2 = lambda wildcards:
            get_sample(
                'fq2_polya_5p',
                search_id='index',
                search_value=wildcards.sample)

    singularity:
        "docker://zavolab/cutadapt:1.16-slim"

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_polya_cutadapt.pe.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "remove_polya_cutadapt.pe.stdout.log")

    shell:
        "(cutadapt \
        -j {threads} \
        --pair-filter=any \
        -m 10 \
        -n 1 \
        -e 0.1 \
        -O 1 \
        -a {params.polya_3_mate1} \
        -g {params.polya_5_mate1} \
        -A {params.polya_3_mate2} \
        -G {params.polya_5_mate2} \
        -o {output.reads1} \
        -p {output.reads2} \
        {input.reads1} \
        {input.reads2}) \
        1> {log.stdout} 2>{log.stderr}"


rule pe_map_genome_star:
    '''
        Map to genome using STAR
    '''
    input:
        index = lambda wildcards:
            os.path.join(
                config["star_indexes"],
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                get_sample(
                    'index_size',
                    search_id='index',
                    search_value=wildcards.sample),
                "STAR_index",
                "chrNameLength.txt"),
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate2.fastq.gz")

    output:
        bam = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.pe.Aligned.sortedByCoord.out.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.pe.Log.final.out")

    params:
        sample_id = "{sample}",
        index = lambda wildcards:
            os.path.join(
                config["star_indexes"],
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                get_sample(
                    'index_size',
                    search_id='index',
                    search_value=wildcards.sample),
                "STAR_index"),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.pe."),
        multimappers = lambda wildcards:
            get_sample(
                'multimappers',
                search_id='index',
                search_value=wildcards.sample),
        soft_clip = lambda wildcards:
            get_sample(
                'soft_clip',
                search_id='index',
                search_value=wildcards.sample),
        pass_mode = lambda wildcards:
            get_sample(
                'pass_mode',
                search_id='index',
                search_value=wildcards.sample),

    singularity:
        "docker://zavolab/star:2.7.3a-slim"

    threads: 12

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "map_genome_star.pe.stderr.log")

    shell:
        "(STAR \
        --runMode alignReads \
        --twopassMode {params.pass_mode} \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.reads1} {input.reads2} \
        --readFilesCommand zcat \
        --outSAMunmapped None  \
        --outFilterMultimapNmax {params.multimappers} \
        --outFilterMultimapScoreRange 0 \
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


rule pe_quantification_salmon:
    '''
        Quantification at transcript and gene level using Salmon
    '''
    input:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate2.fastq.gz"),
        gtf = lambda wildcards:
            get_sample(
                'gtf',
                search_id='index',
                search_value=wildcards.sample),
        index = lambda wildcards:
            os.path.join(
                config["salmon_indexes"],
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                get_sample(
                    'kmer',
                    search_id='index',
                    search_value=wildcards.sample),
                "salmon.idx")

    output:
        gn_estimates = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe",
            "quant.genes.sf"),
        tr_estimates = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe",
            "quant.sf")

    params:
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe"),
        libType = lambda wildcards:
            get_sample(
                'libtype',
                search_id='index',
                search_value=wildcards.sample)

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "genome_quantification_salmon.pe.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "genome_quantification_salmon.pe.stdout.log"),

    threads: 6

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
        -1 {input.reads1} \
        -2 {input.reads2} \
        -o {params.output_dir}; \
        ) 1> {log.stdout} 2> {log.stderr}"


rule pe_genome_quantification_kallisto:
    '''
        Quantification at transcript and gene level using Kallisto
    '''
    input:
        reads1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate1.fastq.gz"),
        reads2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate2.fastq.gz"),
        index = lambda wildcards:
            os.path.join(
                config["kallisto_indexes"],
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                "kallisto.idx")

    output:
        pseudoalignment = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto",
            "{sample}.pe.kallisto.pseudo.sam")

    params:
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto"),
        directionality = lambda wildcards:
            get_sample(
                'kallisto_directionality',
                search_id='index',
                search_value=wildcards.sample)

    singularity:
        "docker://zavolab/kallisto:0.46.1-slim"

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "genome_quantification_kallisto.pe.stderr.log")

    shell:
        "(kallisto quant \
        -i {input.index} \
        -o {params.output_dir} \
        --pseudobam \
        {params.directionality}-stranded \
        {input.reads1} {input.reads2} > {output.pseudoalignment}) \
        2> {log.stderr}"

