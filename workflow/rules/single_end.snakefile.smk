current_rule = 'remove_adapters_cutadapt'
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
        reads = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_adapters_mate1.fastq.gz"))

    params:
        cluster_log_path = config["cluster_log_dir"],
        adapters_3 = lambda wildcards:
            get_sample(
                'fq1_3p',
                search_id='index',
                search_value=wildcards.sample),
        adapters_5 = lambda wildcards:
            get_sample(
                'fq1_5p',
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
            current_rule + ".se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".se.stdout.log")
    shell:
        "(cutadapt \
        -j {threads} \
        -a {params.adapters_3} \
        -g {params.adapters_5} \
        -m 1 \
        {params.additional_params} \
        -o {output.reads} \
        {input.reads}) \
        1> {log.stdout} 2> {log.stderr}"


current_rule = 'remove_polya_cutadapt'
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
        reads = temp(os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.se.remove_polya_mate1.fastq.gz"))

    params:
        cluster_log_path = config["cluster_log_dir"],
        polya_3 = lambda wildcards:
            get_sample(
                'fq1_polya_3p',
                search_id='index',
                search_value=wildcards.sample),
        polya_5 = lambda wildcards:
            get_sample(
                'fq1_polya_5p',
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
            current_rule + ".se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".se.stdout.log")

    shell:
        "(cutadapt \
        -j {threads} \
        -a {params.polya_3} \
        -g {params.polya_5} \
        -m 1 \
        {params.additional_params} \
        -o {output.reads} \
        {input.reads};) \
        1> {log.stdout} 2> {log.stderr}"


current_rule = 'map_genome_star'
rule map_genome_star:
    '''
        Map to genome using STAR
    '''
    input:
        index = lambda wildcards:
            os.path.join(
                config["star_indexes"],
                get_sample('organism', search_id='index', search_value=wildcards.sample),
                get_sample('index_size', search_id='index', search_value=wildcards.sample),
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
            "{sample}.se.Aligned.out.bam"),
        logfile = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.se.Log.final.out")

    shadow: "minimal"

    params:
        cluster_log_path = config["cluster_log_dir"],
        sample_id = "{sample}",
        index = lambda wildcards:
            os.path.abspath(os.path.join(
                config["star_indexes"],
                get_sample('organism', search_id='index', search_value=wildcards.sample),
                get_sample('index_size', search_id='index', search_value=wildcards.sample),
                "STAR_index")),
        outFileNamePrefix = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.se."),
        additional_params = parse_rule_config(
            rule_config,
            current_rule=current_rule,
            immutable=(
                '--genomeDir',
                '--readFilesIn',
                '--readFilesCommand',
                '--outFileNamePrefix',
                '--outSAMattributes',
                '--outStd',
                '--outSAMtype',
                '--outSAMattrRGline',
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
            current_rule + ".se.stderr.log")

    shell:
        "(STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input.reads} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --outSAMattributes All \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outSAMattrRGline ID:rnaseq_pipeline SM:{params.sample_id} \
        {params.additional_params} \
        > {output.bam};) \
        2> {log.stderr}"


current_rule = 'quantification_salmon'
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
                get_sample(
                    'organism',
                    search_id='index',
                    search_value=wildcards.sample),
                get_sample(
                    'kmer',
                    search_id='index',
                    search_value=wildcards.sample),
                "salmon.idx"),
        gtf = lambda wildcards:
            os.path.abspath(get_sample(
                'gtf',
                search_id='index',
                search_value=wildcards.sample))

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
            "quant.sf"),
        meta_info = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se",
            "aux_info",
            "meta_info.json"),
        flenDist = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se",
            "libParams",
            "flenDist.txt")

    shadow: "minimal"

    params:
        cluster_log_path = config["cluster_log_dir"],
        output_dir = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "{sample}.salmon.se"),
        libType = lambda wildcards:
            get_sample(
                'libtype',
                search_id='index',
                search_value=wildcards.sample),
        fraglen = lambda wildcards:
            get_sample(
                'mean',
                search_id='index',
                search_value=wildcards.sample),
        fragsd = lambda wildcards:
            get_sample(
                'sd',
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
                '--unmatedReads',
                '-o',
                )
            )
    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".se.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".se.stdout.log")

    threads: 12

    singularity:
        "docker://quay.io/biocontainers/salmon:1.4.0--h84f40af_1"

    conda:
        os.path.join(workflow.basedir, "envs", "salmon.yaml")

    shell:
        "(salmon quant \
        --libType {params.libType} \
        --threads {threads} \
        --fldMean {params.fraglen} \
        --fldSD {params.fragsd} \
        {params.additional_params} \
        --index {input.index} \
        --geneMap {input.gtf} \
        --unmatedReads {input.reads} \
        -o {params.output_dir};) \
        1> {log.stdout} 2> {log.stderr}"


current_rule = 'genome_quantification_kallisto'
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
            "{sample}.se.kallisto.pseudo.sam"),
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
        fraglen = lambda wildcards:
            get_sample(
                'mean',
                search_id='index',
                search_value=wildcards.sample),
        fragsd = lambda wildcards:
            get_sample(
                'sd',
                search_id='index',
                search_value=wildcards.sample),
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

    threads: 8

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            current_rule + ".se.stderr.log")

    singularity:
        "docker://quay.io/biocontainers/kallisto:0.46.2--h60f4f9f_2"

    conda:
        os.path.join(workflow.basedir, "envs", "kallisto.yaml")

    shell:
        "(kallisto quant \
        --single \
        -i {input.index} \
        -o {params.output_dir} \
        -l {params.fraglen} \
        -s {params.fragsd} \
        -t {threads} \
        {params.directionality} \
        {params.additional_params} \
        --pseudobam \
        {input.reads} > {output.pseudoalignment};) \
        2> {log.stderr}"

