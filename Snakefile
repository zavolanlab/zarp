"""General purpose RNA-Seq analysis pipeline developed by the Zavolan Lab"""
import os
import pandas as pd
import shutil

# Get sample table
samples_table = pd.read_csv(
    config['samples'],
    header=0,
    index_col=0,
    comment='#',
    engine='python',
    sep="\t",
)

# Global config
localrules: start, finish, rename_star_rpm_for_alfa, prepare_multiqc_config


if cluster_config:
    os.makedirs(
        os.path.join(
            os.getcwd(),
            os.path.dirname(cluster_config['__default__']['out']),
        ),
        exist_ok=True)

# Include subworkflows
include: os.path.join("workflow", "rules", "paired_end.snakefile.smk")
include: os.path.join("workflow", "rules", "single_end.snakefile.smk")

rule finish:
    """
        Rule for collecting outputs
    """
    input:
        multiqc_report = os.path.join(
            config['output_dir'],
            "multiqc_summary"),
        bigWig = expand(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "bigWig",
                "{unique_type}",
                "{sample}_{unique_type}_{strand}.bw"),
            sample=samples_table.index.values,
            strand=["plus", "minus"],
            unique_type=["Unique", "UniqueMultiple"]),

        salmon_merge_genes = expand(
            os.path.join(
                config["output_dir"],
                "summary_salmon",
                "quantmerge",
                "genes_{salmon_merge_on}.tsv"),
            salmon_merge_on=["tpm", "numreads"]),

        salmon_merge_transcripts = expand(
            os.path.join(
                config["output_dir"],
                "summary_salmon",
                "quantmerge",
                "transcripts_{salmon_merge_on}.tsv"),
            salmon_merge_on=["tpm", "numreads"]),


rule start:
    '''
       Get samples
    '''
    input:
        reads = lambda wildcards:
            samples_table.loc[wildcards.sample, wildcards.mate],

    output:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "start",
            "{sample}.{mate}.fastq.gz")

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "start_{sample}.{mate}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "start_{sample}.{mate}.stdout.log")

    singularity:
        "docker://bash:5.0.16"

    shell:
        "(cp {input.reads} {output.reads}) \
        1> {log.stdout} 2> {log.stderr} "


rule fastqc:
    '''
        A quality control tool for high throughput sequence data
    '''
    input:
        reads = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "start",
            "{sample}.{mate}.fastq.gz")

    output:
        outdir = directory(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "fastqc",
                "{mate}"))

    threads: 2

    singularity:
        "docker://zavolab/fastqc:0.11.9-slim"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "fastqc_{mate}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "fastqc_{mate}.stdout.log")

    shell:
        "(mkdir -p {output.outdir}; \
        fastqc --outdir {output.outdir} {input.reads}) \
        1> {log.stdout} 2> {log.stderr}"


rule create_index_star:
    """
        Create index for STAR alignments
    """
    input:
        genome = lambda wildcards:
            samples_table['genome']
            [samples_table['organism'] == wildcards.organism]
            [0],
        gtf = lambda wildcards:
            samples_table['gtf']
            [samples_table['organism'] == wildcards.organism]
            [0]

    output:
        chromosome_info = os.path.join(
            config['star_indexes'],
            "{organism}",
            "{index_size}",
            "STAR_index",
            "chrNameLength.txt"),
        chromosomes_names = os.path.join(
            config['star_indexes'],
            "{organism}",
            "{index_size}",
            "STAR_index",
            "chrName.txt")

    params:
        output_dir = os.path.join(
            config['star_indexes'],
            "{organism}",
            "{index_size}",
            "STAR_index"),
        outFileNamePrefix = os.path.join(
            config['star_indexes'],
            "{organism}",
            "{index_size}",
            "STAR_index/STAR_"),
        sjdbOverhang = "{index_size}"

    singularity:
        "docker://zavolab/star:2.7.3a-slim"

    threads: 12

    log:
        stderr = os.path.join(
            config['log_dir'],
            "{organism}_{index_size}_create_index_star.stderr.log"),
        stdout = os.path.join(
            config['log_dir'],
            "{organism}_{index_size}_create_index_star.stdout.log")

    shell:
        "(mkdir -p {params.output_dir}; \
        chmod -R 777 {params.output_dir}; \
        STAR \
        --runMode genomeGenerate \
        --sjdbOverhang {params.sjdbOverhang} \
        --genomeDir {params.output_dir} \
        --genomeFastaFiles {input.genome} \
        --runThreadN {threads} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --sjdbGTFfile {input.gtf}) \
        1> {log.stdout} 2> {log.stderr}"


rule extract_transcriptome:
    """
        Create transcriptome from genome and gene annotations
    """
    input:
        genome = lambda wildcards:
            samples_table['genome'][
                samples_table['organism'] == wildcards.organism][0],
        gtf = lambda wildcards:
            samples_table['gtf'][
                samples_table['organism'] == wildcards.organism][0]

    output:
        transcriptome = os.path.join(
            config['output_dir'],
            "transcriptome",
            "{organism}",
            "transcriptome.fa")

    log:
        stderr = os.path.join(
            config['log_dir'],
            "{organism}_extract_transcriptome.log"),
        stdout = os.path.join(
            config['log_dir'],
            "{organism}_extract_transcriptome.log")

    singularity:
        "docker://zavolab/gffread:0.11.7-slim"

    shell:
        "(gffread \
        -w {output.transcriptome} \
        -g {input.genome} {input.gtf}) \
        1> {log.stdout} 2> {log.stderr}"


rule concatenate_transcriptome_and_genome:
    """
        Concatenate genome and transcriptome
    """
    input:
        transcriptome = os.path.join(
            config['output_dir'],
            "transcriptome",
            "{organism}",
            "transcriptome.fa"),

        genome = lambda wildcards:
            samples_table['genome']
            [samples_table['organism'] == wildcards.organism]
            [0]

    output:
        genome_transcriptome = os.path.join(
            config['output_dir'],
            "transcriptome",
            "{organism}",
            "genome_transcriptome.fa")

    singularity:
        "docker://bash:5.0.16"

    log:
        stderr = os.path.join(
            config['log_dir'],
            "{organism}_concatenate_transcriptome_and_genome.stderr.log")

    shell:
        "(cat {input.transcriptome} {input.genome} \
        1> {output.genome_transcriptome}) \
        2> {log.stderr}"


rule create_index_salmon:
    """
        Create index for Salmon quantification
    """
    input:
        genome_transcriptome = os.path.join(
            config['output_dir'],
            "transcriptome",
            "{organism}",
            "genome_transcriptome.fa"),
        chr_names = lambda wildcards:
            os.path.join(
                config['star_indexes'],
                samples_table["organism"][0],
                str(samples_table["index_size"][0]),
                "STAR_index",
                "chrName.txt")

    output:
        index = directory(
            os.path.join(
                config['salmon_indexes'],
                "{organism}",
                "{kmer}",
                "salmon.idx"))

    params:
        kmerLen = "{kmer}"

    singularity:
        "docker://zavolab/salmon:1.1.0-slim"

    log:
        stderr = os.path.join(
            config['log_dir'],
            "{organism}_{kmer}_create_index_salmon.stderr.log"),
        stdout = os.path.join(
            config['log_dir'],
            "{organism}_{kmer}_create_index_salmon.stdout.log")

    threads: 8

    shell:
        "(salmon index \
        --transcripts {input.genome_transcriptome} \
        --decoys {input.chr_names} \
        --index {output.index} \
        --kmerLen {params.kmerLen} \
        --threads {threads}) \
        1> {log.stdout} 2> {log.stderr}"


rule create_index_kallisto:
    """
        Create index for Kallisto quantification
    """
    input:
        transcriptome = os.path.join(
            config['output_dir'],
            "transcriptome",
            "{organism}",
            "transcriptome.fa")

    output:
        index = os.path.join(
            config['kallisto_indexes'],
            "{organism}",
            "kallisto.idx")

    params:
        output_dir = os.path.join(
            config['kallisto_indexes'],
            "{organism}")

    singularity:
        "docker://zavolab/kallisto:0.46.1-slim"

    log:
        stderr = os.path.join(
            config['log_dir'],
            "{organism}_create_index_kallisto.stderr.log"),
        stdout = os.path.join(
            config['log_dir'],
            "{organism}_create_index_kallisto.stdout.log")

    shell:
        "(mkdir -p {params.output_dir}; \
        chmod -R 777 {params.output_dir}; \
        kallisto index -i {output.index} {input.transcriptome}) \
        1> {log.stdout}  2> {log.stderr}"


rule extract_transcripts_as_bed12:
    """
        Convert transcripts to BED12 format
    """
    input:
        gtf = lambda wildcards:
            samples_table['gtf'][0]

    output:
        bed12 = os.path.join(
            config['output_dir'],
            "full_transcripts_protein_coding.bed")

    singularity:
        "docker://zavolab/zgtf:0.1"

    threads: 1

    log:
        stdout = os.path.join(
            config['log_dir'],
            "extract_transcripts_as_bed12.stdout.log"),
        stderr = os.path.join(
            config['log_dir'],
            "extract_transcripts_as_bed12.stderr.log")

    shell:
        "(gtf2bed12 \
        --gtf {input.gtf} \
        --transcript_type protein_coding \
        --bed12 {output.bed12}); \
        1> {log.stdout} 2> {log.stderr}"


rule index_genomic_alignment_samtools:
    '''
        Index genome bamfile using samtools
    '''
    input:
        bam = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam"),
    output:
        bai = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "map_genome",
            "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam.bai")

    singularity:
        "docker://zavolab/samtools:1.10-slim"

    threads: 1

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "index_genomic_alignment_samtools.{seqmode}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "index_genomic_alignment_samtools.{seqmode}.stdout.log")

    shell:
        "(samtools index {input.bam} {output.bai};) \
        1> {log.stdout} 2> {log.stderr}"


rule calculate_TIN_scores:
    """
        Calculate transcript integrity (TIN) score
    """
    input:
        bam = lambda wildcards:
            expand(
                os.path.join(
                    config['output_dir'],
                    "samples",
                    "{sample}",
                    "map_genome",
                    "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam"),
                sample=wildcards.sample,
                seqmode=samples_table.loc[wildcards.sample, 'seqmode']),
        bai = lambda wildcards:
            expand(
                os.path.join(
                    config['output_dir'],
                    "samples",
                    "{sample}",
                    "map_genome",
                    "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam.bai"),
                sample=wildcards.sample,
                seqmode=samples_table.loc[wildcards.sample, 'seqmode']),
        transcripts_bed12 = os.path.join(
            config['output_dir'],
            "full_transcripts_protein_coding.bed")

    output:
        TIN_score = os.path.join(
            config['output_dir'],
            "samples",
            "{sample}",
            "TIN",
            "TIN_score.tsv")

    params:
        sample = "{sample}"

    log:
        stderr = os.path.join(
            config['log_dir'],
            "samples",
            "{sample}",
            "calculate_TIN_scores.log")

    threads: 8

    singularity:
        "docker://zavolab/tin_score_calculation:0.2.0-slim"

    shell:
        "(tin_score_calculation.py \
        -i {input.bam} \
        -r {input.transcripts_bed12} \
        -c 0 \
        --names {params.sample} \
        -n 100 > {output.TIN_score};) 2> {log.stderr}"


rule merge_TIN_scores:
    """
        Merge TIN scores tables
    """
    input:
        TIN_score = expand(
            os.path.join(
                config['output_dir'],
                "samples",
                "{sample}",
                "TIN",
                "TIN_score.tsv"),
            sample=samples_table.index.values),

    output:
        TIN_scores_merged = os.path.join(
            config['output_dir'],
            "TIN_scores_merged.tsv")

    log:
        stderr = os.path.join(
            config['log_dir'],
            "merge_TIN_scores.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "merge_TIN_scores.stdout.log")

    params:
        TIN_score_merged_paths = ",".join(expand(
            os.path.join(
                config['output_dir'],
                "samples",
                "{sample}",
                "TIN",
                "TIN_score.tsv"),
            zip,
            sample=[i for i in list(samples_table.index.values)],
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)]))

    threads: 1

    singularity:
        "docker://zavolab/tin_score_calculation:0.2.0-slim"

    shell:
        "(tin_score_merge.py \
        --input-files {params.TIN_score_merged_paths} \
        --output-file {output.TIN_scores_merged}) \
        1> {log.stdout} 2> {log.stderr}"


rule plot_TIN_scores:
    """
        Generate TIN scores boxplots
    """
    input:
        TIN_scores_merged = os.path.join(
            config['output_dir'],
            "TIN_scores_merged.tsv"),

    output:
        TIN_boxplot_PNG = os.path.join(
            config['output_dir'],
            "TIN_scores_boxplot_mqc.png"),
        TIN_boxplot_PDF = os.path.join(
            config['output_dir'],
            "TIN_scores_boxplot_mqc.pdf")

    params:
        TIN_boxplot_prefix = os.path.join(
            config['output_dir'],
            "TIN_scores_boxplot_mqc")

    log:
        stderr = os.path.join(
            config['log_dir'],
            "plot_TIN_scores.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "plot_TIN_scores.stdout.log")

    threads: 1

    singularity:
        "docker://zavolab/tin_score_calculation:0.2.0-slim"

    shell:
        "(tin_score_plot.py \
        --input-file {input.TIN_scores_merged} \
        --output-file-prefix {params.TIN_boxplot_prefix}) \
        1> {log.stdout} 2> {log.stderr}"


rule salmon_quantmerge_genes:
    '''
        Merge gene quantifications into a single file
    '''
    input:
        salmon_in = expand(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "{sample}.salmon.{seqmode}",
                "quant.sf"),
            zip,
            sample=samples_table.index.values,
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)])

    output:
        salmon_out = os.path.join(
            config["output_dir"],
            "summary_salmon",
            "quantmerge",
            "genes_{salmon_merge_on}.tsv")

    params:
        salmon_in = expand(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "{sample}.salmon.{seqmode}"),
            zip,
            sample=[i for i in list(samples_table.index.values)],
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)]),
        sample_name_list = expand(
            "{sample}",
            sample=list(samples_table.index.values)),
        salmon_merge_on = "{salmon_merge_on}"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "salmon_quantmerge_genes_{salmon_merge_on}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "salmon_quantmerge_genes_{salmon_merge_on}.stdout.log")

    threads: 1

    singularity:
        "docker://zavolab/salmon:1.1.0-slim"

    shell:
        "(salmon quantmerge \
        --quants {params.salmon_in} \
        --genes \
        --names {params.sample_name_list} \
        --column {params.salmon_merge_on} \
        --output {output.salmon_out};) \
        1> {log.stdout} 2> {log.stderr}"


rule salmon_quantmerge_transcripts:
    '''
        Merge gene quantifications into a single file
    '''
    input:
        salmon_in = expand(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "{sample}.salmon.{seqmode}",
                "quant.sf"),
            zip,
            sample=[i for i in list(samples_table.index.values)],
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)])

    output:
        salmon_out = os.path.join(
            config["output_dir"],
            "summary_salmon",
            "quantmerge",
            "transcripts_{salmon_merge_on}.tsv")

    params:
        salmon_in = expand(
            os.path.join(
                config["output_dir"],
                "samples",
                "{sample}",
                "{sample}.salmon.{seqmode}"),
            zip,
            sample=[i for i in list(samples_table.index.values)],
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)]),

        sample_name_list = expand(
            "{sample}",
            sample=list(samples_table.index.values)),
        salmon_merge_on = "{salmon_merge_on}"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "salmon_quantmerge_transcripts_{salmon_merge_on}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "salmon_quantmerge_transcripts_{salmon_merge_on}.stdout.log")

    threads: 1

    singularity:
        "docker://zavolab/salmon:1.1.0-slim"

    shell:
        "(salmon quantmerge \
        --quants {params.salmon_in} \
        --names {params.sample_name_list} \
        --column {params.salmon_merge_on} \
        --output {output.salmon_out}) \
        1> {log.stdout} 2> {log.stderr}"


rule star_rpm:
    '''
        Create stranded bedgraph coverage with STARs RPM normalisation
    '''
    input:
        bam = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "samples",
                    "{sample}",
                    "map_genome",
                    "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam"),
                sample=wildcards.sample,
                seqmode=samples_table.loc[wildcards.sample, 'seqmode']),
        bai = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "samples",
                    "{sample}",
                    "map_genome",
                    "{sample}.{seqmode}.Aligned.sortedByCoord.out.bam.bai"),
                sample=wildcards.sample,
                seqmode=samples_table.loc[wildcards.sample, 'seqmode']),

    output:
        str1 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "STAR_coverage",
            "{sample}_Signal.Unique.str1.out.bg"),
        str2 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "STAR_coverage",
            "{sample}_Signal.UniqueMultiple.str1.out.bg"),
        str3 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "STAR_coverage",
            "{sample}_Signal.Unique.str2.out.bg"),
        str4 = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "STAR_coverage",
            "{sample}_Signal.UniqueMultiple.str2.out.bg")

    params:
        out_dir = lambda wildcards, output:
            os.path.dirname(output.str1),
        prefix = lambda wildcards, output:
            os.path.join(
                os.path.dirname(output.str1),
                str(wildcards.sample) + "_"),
        stranded = "Stranded"

    singularity:
        "docker://zavolab/star:2.7.3a-slim"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "star_rpm.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "star_rpm.stdout.log")

    threads: 4

    shell:
        "(mkdir -p {params.out_dir}; \
        chmod -R 777 {params.out_dir}; \
        STAR \
        --runMode inputAlignmentsFromBAM \
        --runThreadN {threads} \
        --inputBAMfile {input.bam} \
        --outWigType bedGraph \
        --outWigStrand {params.stranded} \
        --outWigNorm RPM \
        --outFileNamePrefix {params.prefix}) \
        1> {log.stdout} 2> {log.stderr}"


rule rename_star_rpm_for_alfa:
    input:
        plus = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "samples",
                    "{sample}",
                    "STAR_coverage",
                    "{sample}_Signal.{unique}.{plus}.out.bg"),
                sample=wildcards.sample,
                unique=wildcards.unique,
                plus=samples_table.loc[wildcards.sample, 'alfa_plus']),

        minus = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "samples",
                    "{sample}",
                    "STAR_coverage",
                    "{sample}_Signal.{unique}.{minus}.out.bg"),
                sample=wildcards.sample,
                unique=wildcards.unique,
                minus=samples_table.loc[wildcards.sample, 'alfa_minus'])

    output:
        plus = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.{unique}.plus.bg"),
        minus = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.{unique}.minus.bg")

    params:
        orientation = lambda wildcards:
            samples_table.loc[wildcards.sample, "kallisto_directionality"]

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "rename_star_rpm_for_alfa__{unique}.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "rename_star_rpm_for_alfa__{unique}.stdout.log")

    singularity:
        "docker://bash:5.0.16"

    shell:
        "(cp {input.plus} {output.plus}; \
         cp {input.minus} {output.minus};) \
         1>{log.stdout} 2>{log.stderr}"


rule generate_alfa_index:
    ''' Generate ALFA index files from sorted GTF file '''
    input:
        gtf = lambda wildcards:
            samples_table["gtf"]
            [samples_table["organism"] == wildcards.organism][0],
        chr_len = os.path.join(
            config["star_indexes"],
            "{organism}",
            "{index_size}",
            "STAR_index",
            "chrNameLength.txt"),

    output:
        index_stranded = os.path.join(
            config["alfa_indexes"],
            "{organism}",
            "{index_size}",
            "ALFA",
            "sorted_genes.stranded.ALFA_index"),
        index_unstranded = os.path.join(
            config["alfa_indexes"],
            "{organism}",
            "{index_size}",
            "ALFA",
            "sorted_genes.unstranded.ALFA_index")

    params:
        genome_index = "sorted_genes",
        out_dir = lambda wildcards, output:
            os.path.dirname(output.index_stranded)

    threads: 4

    singularity:
        "docker://zavolab/alfa:1.1.1-slim"

    log:
        os.path.join(
            config["log_dir"],
            "{organism}_{index_size}_generate_alfa_index.log")

    shell:
        "(alfa -a {input.gtf} \
        -g {params.genome_index} \
        --chr_len {input.chr_len} \
        -p {threads} \
        -o {params.out_dir}) &> {log}"


rule alfa_qc:
    '''
        Run ALFA from stranded bedgraph files
    '''
    input:
        plus = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.{unique}.plus.bg"),
        minus = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.{unique}.minus.bg"),
        gtf = lambda wildcards:
            os.path.join(
                config["alfa_indexes"],
                samples_table.loc[wildcards.sample, "organism"],
                str(samples_table.loc[wildcards.sample, "index_size"]),
                "ALFA",
                "sorted_genes.stranded.ALFA_index")

    output:
        biotypes = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "ALFA_plots.Biotypes.pdf"),
        categories = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "ALFA_plots.Categories.pdf"),
        table = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.ALFA_feature_counts.tsv")

    params:
        out_dir = lambda wildcards, output:
            os.path.dirname(output.biotypes),
        plus = lambda wildcards, input:
            os.path.basename(input.plus),
        minus = lambda wildcards, input:
            os.path.basename(input.minus),
        alfa_orientation = lambda wildcards:
            [samples_table.loc[
                wildcards.sample, "alfa_directionality"]],
        genome_index = lambda wildcards, input:
            os.path.abspath(
                os.path.join(
                    os.path.dirname(input.gtf),
                    "sorted_genes")),
        name = "{sample}"

    singularity:
        "docker://zavolab/alfa:1.1.1-slim"

    log:
        os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "alfa_qc.{unique}.log")

    shell:
        "(cd {params.out_dir}; \
        alfa \
        -g {params.genome_index} \
        --bedgraph {params.plus} {params.minus} {params.name} \
        -s {params.alfa_orientation}) &> {log}"


# cd {params.out_dir};
rule alfa_qc_all_samples:
    '''
        Run ALFA from stranded bedgraph files on all samples
    '''
    input:
        tables = lambda wildcards:
            expand(
                os.path.join(
                    config["output_dir"],
                    "samples",
                    "{sample}",
                    "ALFA",
                    "{unique}",
                    "{sample}.ALFA_feature_counts.tsv"),
                sample=samples_table.index.values,
                unique=wildcards.unique)
    output:
        biotypes = os.path.join(
            config["output_dir"],
            "ALFA",
            "{unique}",
            "ALFA_plots.Biotypes.pdf"),
        categories = os.path.join(
            config["output_dir"],
            "ALFA",
            "{unique}",
            "ALFA_plots.Categories.pdf")

    params:
        out_dir = lambda wildcards, output:
            os.path.dirname(output.biotypes)

    log:
        os.path.join(
            config["log_dir"],
            "alfa_qc_all_samples.{unique}.log")

    singularity:
        "docker://zavolab/alfa:1.1.1-slim"

    shell:
        "(alfa -c {input.tables} -o {params.out_dir}) &> {log}"


rule alfa_concat_results:
    input:
        expand(
            os.path.join(
                config["output_dir"],
                "ALFA",
                "{unique}",
                "ALFA_plots.{annotation}.pdf"),
            unique=["Unique", "UniqueMultiple"],
            annotation=["Categories", "Biotypes"])

    output:
        os.path.join(
            config["output_dir"],
            "ALFA",
            "ALFA_plots_mqc.png")

    params:
        density = 300

    log:
        os.path.join(
            config["log_dir"],
            "alfa_qc_all_samples.concat.log")

    singularity:
        "docker://zavolab/imagemagick:7.0.8"

    shell:
        "(convert -append -density {params.density} \
            {input} {output}) &> {log}"


rule prepare_multiqc_config:
    '''
        Prepare config for the MultiQC
    '''
    input:
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "scripts",
            "zarp_multiqc_config.py")

    output:
        multiqc_config = os.path.join(
            config["output_dir"],
            "multiqc_config.yaml")

    params:
        logo_path = config['report_logo'],
        multiqc_intro_text = config['report_description'],
        url = config['report_url']

    log:
        stderr = os.path.join(
            config["log_dir"],
            "prepare_multiqc_config.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "prepare_multiqc_config.stdout.log")

    shell:
        "(python {input.script} \
        --config {output.multiqc_config} \
        --intro-text '{params.multiqc_intro_text}' \
        --custom-logo {params.logo_path} \
        --url '{params.url}') \
        1> {log.stdout} 2> {log.stderr}"


rule multiqc_report:
    '''
        Create report with MultiQC
    '''
    input:
        fastqc_se = expand(
            os.path.join(
                config['output_dir'],
                "samples",
                "{sample}",
                "fastqc",
                "{mate}"),
            sample=samples_table.index.values,
            mate="fq1"),

        fastqc_pe = expand(
            os.path.join(
                config['output_dir'],
                "samples",
                "{sample}",
                "fastqc",
                "{mate}"),
            sample=[i for i in list(
                samples_table[samples_table['seqmode'] == 'pe'].index.values)],
            mate="fq2"),

        pseudoalignment = expand(
            os.path.join(
                config['output_dir'],
                "samples",
                "{sample}",
                "quant_kallisto",
                "{sample}.{seqmode}.kallisto.pseudo.sam"),
            zip,
            sample=[i for i in list(samples_table.index.values)],
            seqmode=[samples_table.loc[i, 'seqmode']
                     for i in list(samples_table.index.values)]),

        TIN_boxplot_PNG = os.path.join(
            config['output_dir'],
            "TIN_scores_boxplot_mqc.png"),

        TIN_boxplot_PDF = os.path.join(
            config['output_dir'],
            "TIN_scores_boxplot_mqc.pdf"),

        alfa_concat_out = os.path.join(
            config["output_dir"],
            "ALFA",
            "ALFA_plots_mqc.png"),

        multiqc_config = os.path.join(
            config["output_dir"],
            "multiqc_config.yaml")

    output:
        multiqc_report = directory(
            os.path.join(
                config["output_dir"],
                "multiqc_summary"))

    params:
        results_dir = os.path.join(
            config["output_dir"]),
        log_dir = config["log_dir"]

    log:
        stderr = os.path.join(
            config["log_dir"],
            "multiqc_report.stderr.log"),
        stdout = os.path.join(
            config["log_dir"],
            "multiqc_report.stdout.log")

    singularity:
        "docker://ewels/multiqc:1.7"

    shell:
        "(multiqc \
        --outdir {output.multiqc_report} \
        --config {input.multiqc_config} \
        {params.results_dir} \
        {params.log_dir};) \
        1> {log.stdout} 2> {log.stderr}"


rule sort_bed_4_big:
    '''
        sort bedGraphs in order to work with bedGraphtobigWig
    '''
    input:
        bg = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "ALFA",
            "{unique}",
            "{sample}.{unique}.{strand}.bg")

    output:
        sorted_bg = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "bigWig",
            "{unique}",
            "{sample}_{unique}_{strand}.sorted.bg")

    singularity:
        "docker://cjh4zavolab/bedtools:2.27"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "sort_bg_{unique}_{strand}.stderr.log")

    shell:
        "(sortBed \
         -i {input.bg} \
         > {output.sorted_bg};) 2> {log.stderr}"

rule prepare_bigWig:
    '''
        bedGraphtobigWig, for viewing in genome browsers
    '''
    input:
        sorted_bg = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "bigWig",
            "{unique}",
            "{sample}_{unique}_{strand}.sorted.bg"),
        chr_sizes = lambda wildcards:
            os.path.join(
                config['star_indexes'],
                samples_table.loc[wildcards.sample, "organism"],
                str(samples_table.loc[wildcards.sample, "index_size"]),
                "STAR_index",
                "chrNameLength.txt")

    output:
        bigWig = os.path.join(
            config["output_dir"],
            "samples",
            "{sample}",
            "bigWig",
            "{unique}",
            "{sample}_{unique}_{strand}.bw")

    singularity:
        "docker://zavolab/bedgraphtobigwig:4-slim"

    log:
        stderr = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "bigwig_{unique}_{strand}.stderr.log"),

        stdout = os.path.join(
            config["log_dir"],
            "samples",
            "{sample}",
            "bigwig_{unique}_{strand}.stdout.log")

    shell:
        "(bedGraphToBigWig \
         {input.sorted_bg} \
         {input.chr_sizes} \
         {output.bigWig};) \
         1> {log.stdout} 2> {log.stderr}"
