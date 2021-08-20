current_rule = 'pe_remove_adapters_cutadapt'
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
        reads1 = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate1.fastq.gz")),
        reads2 = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_adapters_mate2.fastq.gz"))

    params:
        cluster_log_path = config["cluster_log_dir"],
        adapter_3_mate1 = lambda wildcards:
            get_sample('fq1_3p', search_id='index', search_value=wildcards.sample),
        adapter_5_mate1 = lambda wildcards:
            get_sample('fq1_5p', search_id='index', search_value=wildcards.sample),
        adapter_3_mate2 = lambda wildcards:
            get_sample('fq2_3p', search_id='index', search_value=wildcards.sample),
        adapter_5_mate2 = lambda wildcards:
            get_sample('fq2_5p', search_id='index', search_value=wildcards.sample),
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '-a',
                '-A',
                '-g',
                '-G',
                '-o',
                '-p',
                )
            )

    singularity:
        "docker://quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1"

    conda:
        os.path.join(workflow.basedir, "envs", "cutadapt.yaml")

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stdout.log")

    shell:
        "(cutadapt \
        -j {threads} \
        -a {params.adapter_3_mate1} \
        -g {params.adapter_5_mate1} \
        -A {params.adapter_3_mate2} \
        -G {params.adapter_5_mate2} \
        -m 1 \
        {params.additional_params} \
        -o {output.reads1} \
        -p {output.reads2} \
        {input.reads1} \
        {input.reads2};) \
        1> {log.stdout} 2>{log.stderr}"


current_rule = 'pe_remove_polya_cutadapt'
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
        reads1 = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate1.fastq.gz")),
        reads2 = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.pe.remove_polya_mate2.fastq.gz"))

    params:
        cluster_log_path = config["cluster_log_dir"],
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
                search_value=wildcards.sample),
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '-a',
                '-A',
                '-g',
                '-G',
                '-o',
                '-p',
                )
            )

    singularity:
        "docker://quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1"

    conda:
        os.path.join(workflow.basedir, "envs", "cutadapt.yaml")

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stdout.log")

    shell:
        "(cutadapt \
        -j {threads} \
        -a {params.polya_3_mate1} \
        -g {params.polya_5_mate1} \
        -A {params.polya_3_mate2} \
        -G {params.polya_5_mate2} \
        -m 1 \
        {params.additional_params} \
        -o {output.reads1} \
        -p {output.reads2} \
        {input.reads1} \
        {input.reads2}) \
        1> {log.stdout} 2>{log.stderr}"


current_rule = 'pe_map_genome_star'
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

    shadow: "minimal"

    params:
        cluster_log_path = config["cluster_log_dir"],
        sample_id = "{sample}",
        index = lambda wildcards:
            os.path.abspath(os.path.join(
                config["star_indexes"],
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                get_sample(
                    'index_size',
                    search_id='index',
                    search_value=wildcards.sample),
                    "STAR_index")),
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
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '--twopassMode',
                '--genomeDir',
                '--readFilesIn',
                '--readFilesCommand',
                '--outFilterMultimapNmax',
                '--outFileNamePrefix',
                '--outSAMattributes',
                '--outStd',
                '--outSAMtype',
                '--outSAMattrRGline',
                '--alignEndsType',
                )
            )

    singularity:
        "docker://quay.io/biocontainers/star:2.7.8a--h9ee0642_1"

    conda:
        os.path.join(workflow.basedir, "envs", "STAR.yaml")

    threads: 12

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stderr.log")

    shell:
        "(STAR \
        --twopassMode {params.pass_mode} \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.reads1} {input.reads2} \
        --readFilesCommand zcat \
        --outFilterMultimapNmax {params.multimappers} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_SortedByCoordinate \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:rnaseq_pipeline SM:{params.sample_id} \
        --alignEndsType {params.soft_clip} \
        {params.additional_params} \
        > {output.bam};) \
        2> {log.stderr}"


current_rule = 'pe_quantification_salmon'
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
            os.path.abspath(get_sample(
                'gtf',
                search_id='index',
                search_value=wildcards.sample)),
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
            "quant.sf"),
        meta_info = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe",
            "aux_info",
            "meta_info.json"),
        flenDist = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe",
            "libParams",
            "flenDist.txt")

    shadow: "minimal"

    params:
        cluster_log_path = config["cluster_log_dir"],
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.pe"),
        libType = lambda wildcards:
            get_sample(
                'libtype',
                search_id='index',
                search_value=wildcards.sample),
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '--libType',
                '--fldMean',
                '--fldSD',
                '--index',
                '--geneMap',
                '-1',
                '-2',
                '-o',
                )
            )

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stdout.log"),

    threads: 6

    singularity:
        "docker://quay.io/biocontainers/salmon:1.4.0--h84f40af_1"

    conda:
        os.path.join(workflow.basedir, "envs", "salmon.yaml")

    shell:
        "(salmon quant \
        --libType {params.libType} \
        --threads {threads} \
        {params.additional_params} \
        --index {input.index} \
        --geneMap {input.gtf} \
        -1 {input.reads1} \
        -2 {input.reads2} \
        -o {params.output_dir}; \
        ) 1> {log.stdout} 2> {log.stderr}"


current_rule = 'pe_genome_quantification_kallisto'
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
            "{sample}.pe.kallisto.pseudo.sam"),
        abundances = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto",
            "abundance.h5")

    shadow: "minimal"

    params:
        cluster_log_path = config["cluster_log_dir"],
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "quant_kallisto"),
        directionality = lambda wildcards:
            get_directionality(get_sample(
                    'libtype',
                    search_id='index',
                    search_value=wildcards.sample),"kallisto"),
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '--single',
                '-i',
                '-o',
                '-l',
                '-s',
                '--pseudobam',
                '--fr-stranded',
                '--rf-stranded',
                )
            )


    singularity:
        "docker://quay.io/biocontainers/kallisto:0.46.2--h60f4f9f_2"

    conda:
        os.path.join(workflow.basedir, "envs", "kallisto.yaml")

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".stderr.log")

    shell:
        "(kallisto quant \
        -i {input.index} \
        -o {params.output_dir} \
        -t {threads} \
        {params.directionality} \
        {params.additional_params} \
        --pseudobam \
        {input.reads1} {input.reads2} > {output.pseudoalignment}) \
        2> {log.stderr}"

