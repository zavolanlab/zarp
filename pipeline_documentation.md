# ZARP: workflow documentation

This document describes the individual steps of the workflow. For instructions
on installation and usage please see [here](README.md).

## Table of Contents

- [**Third-party software used**](#third-party-software-used)
- [**Description of workflow steps**](#description-of-workflow-steps)
  - [**Rule graph**](#rule-graph)
  - [**Preparatory**](#preparatory)
    - [**Read samples table**](#read-samples-table)
    - [**Create log directories**](#create-log-directories)
  - [**Sequencing mode-independent**](#sequencing-mode-independent)
    - [**start**](#start)
    - [**create_index_star**](#create_index_star)
    - [**extract_transcriptome**](#extract_transcriptome)
    - [**concatenate_transcriptome_and_genome**](#concatenate_transcriptome_and_genome)
    - [**create_index_salmon**](#create_index_salmon)
    - [**create_index_kallisto**](#create_index_kallisto)
    - [**extract_transcripts_as_bed12**](#extract_transcripts_as_bed12)
    - [**fastqc**](#fastqc)
    - [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
    - [**star_rpm**](#star_rpm)
    - [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
    - [**sort_bed_4_big**](#sort_bed_4_big)
    - [**prepare_bigWig**](#prepare_bigwig)
    - [**calculate_TIN_scores**](#calculate_tin_scores)
    - [**merge_TIN_scores**](#merge_tin_scores)
    - [**plot_TIN_scores**](#plot_tin_scores)
    - [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
    - [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)
    - [**generate_alfa_index**](#generate_alfa_index)
    - [**alfa_qc**](#alfa_qc)
    - [**alfa_qc_all_samples**](#alfa_qc_all_samples)
    - [**alfa_concat_results**](#alfa_concat_results)
    - [**prepare_multiqc_config**](#prepare_multiqc_config)
    - [**multiqc_report**](#multiqc_report)
    - [**finish**](#finish)
  - [**Sequencing mode-specific**](#sequencing-mode-specific)
    - [**remove_adapters_cutadapt**](#remove_adapters_cutadapt)
    - [**remove_polya_cutadapt**](#remove_polya_cutadapt)
    - [**map_genome_star**](#map_genome_star)
    - [**quantification_salmon**](#quantification_salmon)
    - [**genome_quantification_kallisto**](#genome_quantification_kallisto)

## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **ALFA** | [MIT][license-mit] | _"**A**nnotation **L**andscape **F**or **A**ligned reads"_ - _"[...] provides a global overview of features distribution composing NGS dataset(s)"_ | [code][code-alfa] / [manual][docs-alfa] / [publication][pub-alfa] |
| **bedGraphToBigWig** | [MIT][license-mit] | _"Convert a bedGraph file to bigWig format"_ | [code][code-bedgraphtobigwig] / [manual][code-bedgraphtobigwig] |
| **bedtools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][code-bedtools] |
| **cutadapt** | [MIT][license-mit] | _"[...] finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads"_ | [code][code-cutadapt] / [manual][docs-cutadapt] / [publication][pub-cutadapt] |
| **gffread** | [MIT][license-mit] | _"[...] validate, filter, convert and perform various other operations on GFF files"_ | [code][code-gffread] / [manual][docs-gffread] |
| **FastQC** | [GPLv3][license-gpl3] | _"A quality control analysis tool for high throughput sequencing data"_ | [code][code-fastqc] / [manual][docs-fastqc] |
| **ImageMagick** | [custom][license-imagemagick]^ | _"[...] create, edit, compose, or convert bitmap images"_ | [code][code-imagemagick] / [manual][docs-imagemagick] |
| **kallisto** | [BSD-2][license-bsd2] | _"[...] program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads"_ | [code][code-kallisto] / [manual][docs-kallisto] / [publication][pub-kallisto] |
| **MultiQC** | [GPLv3][license-gpl3] | _"Aggregate results from bioinformatics analyses across many samples into a single report"_ | [code][code-multiqc] / [manual][docs-multiqc] / [publication][pub-multiqc] |
| **RSeqC** | [GPLv3][license-gpl3] | _"[...] comprehensively evaluate different aspects of RNA-seq experiments, such as sequence quality, GC bias, polymerase chain reaction bias, nucleotide composition bias, sequencing depth, strand specificity, coverage uniformity and read distribution over the genome structure."_ | [code][code-rseqc] / [manual][docs-rseqc] / [publication][pub-rseqc] |
| **Salmon** | [GPLv3][license-gpl3] | _"Highly-accurate & wicked fast transcript-level quantification from RNA-seq reads using selective alignment"_ | [code][code-salmon] / [manual][docs-salmon] / [publication][pub-salmon] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |
| **STAR** | [MIT][license-mit] | _"**S**pliced **T**ranscripts **A**lignment to a **R**eference"_ - _"RNA-seq aligner"_ | [code][code-star] / [manual][docs-star] / [publication][pub-star] |

^ compatible with [GPLv3][license-gpl3]

## Description of workflow steps

> The workflow consists of three Snakemake files: A main `Snakefile` and an
> individual Snakemake file for each sequencing mode (single-end and
> paired-end), as parameters for some tools differ between sequencing modes.
> The main `Snakefile` contains general steps for the creation of indices and
> other required files derived from the annotations, steps that are applicable
> to both sequencing modes, and steps that deal with gathering, summarizing or
> combining results. Individual steps of the workflow are described briefly, and
> links to the respective software manuals are given. Parameters that can be
> modified by the user (via the samples table) are also described. Descriptions
> for steps for which individual "rules" exist for single- and paired-end
> sequencing libraries are combined, and only differences between the modes are
> highlighted.

### Rule graph

![rule_graph][rule-graph]

Visual representation of workflow. Automatically prepared with
[Snakemake][docs-snakemake].

### Preparatory

#### Read sample table

##### Requirements

- tab-separated values (`.tsv`) file
- first row has to contain parameter names as in [`samples.tsv`](tests/input_files/samples.tsv)
- first column used as sample identifiers

Parameter name | Description | Data type(s)
--- | --- | ---
sample | Descriptive sample name | `str`
seqmode | Required for various steps of the workflow. One of `pe` (for paired-end libraries) or `se` (for single-end libraries). | `str`
fq1 | Path of library file in `.fastq.gz` format (or mate 1 read file for paired-end libraries) | `str`
index_size | Required for [STAR](#third-party-software-used). Ideally the maximum read length minus 1 (`max(ReadLength)-1`). Values lower than maximum read length may result in lower mapping accuracy, while higher values may result in longer processing times. | `int`
kmer | Required for [Salmon](#third-party-software-used). Default value of 31 usually works fine for reads of 75 bp or longer. Consider using lower values of poor mapping is observed. | `int`
fq2 | Path of mate 2 read file in `.fastq.gz` format. Value ignored for for single-end libraries. | `str`
fq1_3p | Required for [Cutadapt](#third-party-software-used). 3' adapter of mate 1. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. | `str`
fq1_5p | Required for [Cutadapt](#third-party-software-used). 5' adapter of mate 1. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. | `str`
fq2_3p | Required for [Cutadapt](#third-party-software-used). 3' adapter of mate 2. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. Value ignored for single-end libraries. | `str`
fq2_5p | Required for [Cutadapt](#third-party-software-used). 5' adapter of mate 2. Use value such as `XXXXXXXXXXXXXXX` if no adapter present or if no trimming is desired. Value ignored for single-end libraries. | `str`
organism | Name or identifier of organism or organism-specific genome resource version. Has to correspond to the naming of provided genome and gene annotation files and directories, like "ORGANISM" in the path below. Example: `GRCh38` | `str`
gtf | Required for [STAR](#third-party-software-used). Path to gene annotation `.fa` file. File needs to be in subdirectory corresponding to `organism` field. Example: `/path/to/GRCh38/gene_annotations.gtf` | `str`
gtf_filtered | Required for [Salmon](#third-party-software-used). Path to filtered gene annotation `.gtf` file. File needs to be in subdirectory corresponding to `organism` field. Example: `/path/to/GRCh38/gene_annotations.filtered.gtf` | `str`
genome | Required for [STAR](#third-party-software-used). Path to genome `.fa` file. File needs to be in subdirectory corresponding to `organism` field. Example: `/path/to/GRCh38/genome.fa` | `str`
sd | Required for [kallisto](#third-party-software-used) and [Salmon](#third-party-software-used), but only for single-end libraries. Estimated standard deviation of fragment length distribution. Can be assessed from, e.g., BioAnalyzer profiles. Value ignored for paired-end libraries. | `int`
mean | Required for [kallisto](#third-party-software-used) and [Salmon](#third-party-software-used), but only for single-end libraries. Estimated mean of fragment length distribution. Can be assessed, e.g., from BioAnalyzer profiles. Value ignored for paired-end libraries. | `int`
multimappers | Required for [STAR](#third-party-software-used). Maximum number of multiple alignments allowed for a read; if exceeded, the read is considered unmapped. | `int`
soft_clip | Required for [STAR](#third-party-software-used). One of `Local` (standard local alignment with soft-clipping allowed) or `EndToEnd` (force end-to-end read alignment, do not soft-clip). | `str`
pass_mode | Required for [STAR](#third-party-software-used). One of `None` (1-pass mapping) or `Basic` (basic 2-pass mapping, with all 1st-pass junctions inserted into the genome indices on the fly). | `str`
libtype | Required for [Salmon](#third-party-software-used). See [Salmon manual][docs-salmon] for allowed values. If in doubt, enter `A` to automatically infer the library type. | `str`
kallisto_directionality | Required for [kallisto](#third-party-software-used) and [ALFA](#third-party-software-used). One of `--fr-stranded` (strand-specific reads, first read forward) and `--rf-stranded` (strand-specific reads, first read reverse) | `str`
fq1_polya3p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 3' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. | `str`
fq1_polya5p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 5' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. | `str`
fq2_polya3p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 3' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. Value ignored for single-end libraries. | `str`
fq2_polya5p | Required for [Cutadapt](#third-party-software-used). Stretch of `A`s or `T`s, depending on read orientation. Trimmed from the 5' end of the read. Use value such as `XXXXXXXXXXXXXXX` if no poly(A) stretch present or if no trimming is desired. Value ignored for single-end libraries. | `str`

#### Create log directories

Sets up logging directories for the workflow run environment. Vanilla Python
statement, not a Snakemake rule.

### Sequencing mode-independent

#### `start`

Copy and rename read files.

> Local rule

- **Input**
  - Reads file (`.fastq.gz`)
- **Output**
  - Reads file, copied, renamed (`.fastq.gz`); used in [**fastqc**](#fastqc) and
    [**remove_adapters_cutadapt**](#remove_adapters_cutadapt)

#### `create_index_star`

Create index for [**STAR**](#third-party-software-used) short read aligner.

> Indices need only be generated once for each combination of genome, set of
> annotations and index size.

- **Input**
  - Genome sequence file (`.fasta`)
  - Gene annotation file (`.gtf`)
- **Parameters**
  - `--sjdbOverhang`: maximum read length - 1; lower values may reduce accuracy,
    higher values may increase STAR runtime; specify in sample table column
    `index_size`
- **Output**
  - STAR index; used in [**map_genome_star**](#map_genome_star)
  - Index includes files:
    - Chromosome length table `chrNameLength.txt`; used in
      [**generate_alfa_index**](#generate_alfa_index) and
      [**prepare_bigWig**](#prepare_bigwig)
    - Chromosome name list `chrName.txt`; used in
      [**create_index_salmon**](#create_index_salmon)

#### `extract_transcriptome`

Create transcriptome from genome and gene annotations with
[**gffread**](#third-party-software-used).

- **Input**
  - Genome sequence file (`.fasta`)
  - Gene annotation file (`.gtf`)
- **Output**
  - Transcriptome sequence file (`.fasta`); used in
    [**concatenate_transcriptome_and_genome**](#concatenate_transcriptome_and_genome)
    and [**create_index_kallisto**](#create_index_kallisto)

#### `concatenate_transcriptome_and_genome`

Concatenate reference transcriptome and genome.

> Required by [Salmon][docs-salmon-selective-alignment]

- **Input**
  - Genome sequence file (`.fasta`)
  - Transcriptome sequence file (`.fasta`); from
    [**extract_transcriptome**](#extract_transcriptome)
- **Output**
  - Transcriptome genome reference file (`.fasta`); used in
    [**create_index_salmon**](#create_index_salmon)

#### `create_index_salmon`

Create index for [**Salmon**](#third-party-software-used) quantification.

> Required if Salmon is to be used in "mapping-based" mode.
>create_index_salmon
> Index is built using an auxiliary k-mer hash over k-mers of a specified
> length (default: 31). While the mapping algorithms will make use of
> arbitrarily long matches between the query and reference, the selected k-mer
> size will act as the minimum acceptable length for a valid match. Thus, a
> smaller value of k may slightly improve sensitivty. Empirically, the default
> value of k = 31 seems to work well for reads of 75 bp or longer. For shorter
> reads, consider using a smaller k.

- **Input**
  - Transcriptome genome reference file (`.fasta`); from
    [**concatenate_transcriptome_and_genome**](#concatenate_transcriptome_and_genome)
  - Chromosome name list `chrName.txt`; from
    [**create_index_star**](#create_index_star)
- **Parameters**
  - `--kmerLen`: k-mer length; specify in sample table column `kmer`
- **Output**
  - Salmon index; used in [**quantification_salmon**](#quantification_salmon)

#### `create_index_kallisto`

Create index for [**kallisto**](#third-party-software-used) quantification.

> The default kmer size of 31 is used in this workflow and is not configurable
> by the user.

- **Input**
  - Transcriptome sequence file (`.fasta`); from
    [**extract_transcriptome**](#extract_transcriptome)
- **Output**
  - kallisto index; used in
    [**genome_quantification_kallisto**](#genome_quantification_kallisto)

#### `extract_transcripts_as_bed12`

Convert transcripts from `.gtf` to extended 12-column `.bed` format with
[custom-script][custom-script-gtf-to-bed12].

- **Input**
  - Gene annotation file (`.gtf`)
- **Output**
  - Transcript annotations file (12-column `.bed`); used in
    [**calculate_TIN_scores**](#calculate_tin_scores)

#### `fastqc`

Prepare quality control report for reads library with
[**FastQC**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from [**start**](#start)
- **Output**
  - FastQC output directory with report (`.txt`) and figures (`.png`); used in
    [**multiqc_report**](#multiqc_report)

#### `index_genomic_alignment_samtools`

Index BAM file with [**SAMtools**](#third-party-software-used).

> Indexing a genome sorted BAM file allows one to quickly extract alignments
> overlapping particular genomic regions. Moreover, indexing is required by
> genome viewers such as IGV so that the viewers can quickly display alignments
> in a genomic region of interest.

- **Input**
  - Alignemnts file (`.bam`); from [**map_genome_star**](#map_genome_star)
- **Output**
  - BAM index file (`.bam.bai`); used in [**star_rpm**](#star_rpm) &
    [**calculate_TIN_scores**](#calculate_tin_scores)

#### `star_rpm`

Create stranded bedGraph coverage (`.bg`) with
[**STAR**](#third-party-software-used).

> Uses STAR's [RPM normalization][docs-star-rpm-norm] functionality.
>
> STAR RPM uses SAM flags to correctly tell where the read and its mate mapped
> to. That is, if mate 1 is mapped to the plus strand, then mate 2 is mapped to
> the minus strand, and STAR will count mate 1 and mate 2 to the plus strand.
> This is in contrast to `bedtools genomecov -bg -split`, where a mate is
> assigned to a strand irrespective of its corresponding mate.

- **Input**
  - Alignments file (`.bam`); from [**map_genome_star**](#map_genome_star)
  - BAM index file (`.bam.bai`); from
    [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
- **Output**
  - Coverage file (`.bg`); used in [**multiqc_report**](#multiqc_report) and
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
- **Non-configurable & non-default**
  - `--outWigStrans="Stranded"`
  - `--outWigNorm="RPM"`

#### `rename_star_rpm_for_alfa`

Rename and copy stranded bedGraph coverage tracks such that they comply with
[**ALFA**](#third-party-software-used).

> Local rule
>
> Renaming to `plus.bg` and `minus.bg` depends on library orientation, which is
> provided by user in sample table column `kallisto_directionality`.

- **Input**
  - Coverage file (`.bg`); from [**star_rpm**](#star_rpm)
- **Output**
  - Coverage file, renamed (`.bg`); used in [**alfa_qc**](#alfa_qc) and
    [**sort_bed_4_big**](#sort_bed_4_big)

#### `sort_bed_4_big`

Sort bedGraph files with [**bedtools**](#third-party-software-used).

- **Input**
  - Coverage files, renamed (`.bg`); from
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
- **Output**
  - Coverage files, sorted (`.bg`); used in [**prepare_bigWig**](#prepare_bigwig)

#### `prepare_bigWig`

Generate bigWig from bedGraph with
[**bedGraphToBigWig**](#third-party-software-used).

- **Input**
  - Coverage files, sorted (`.bg`); from [**sort_bed_4_big**](#sort_bed_4_big)
  - Chromosome length table `chrNameLength.txt`; from
    [**create_index_star**](#create_index_star)
- **Output**
  - Coverage files, one per strand and sample (`.bw`); used in
    [**finish**](#finish)

#### `calculate_TIN_scores`

Calculates the Transcript Integrity Number (TIN) for each transcript with
[custom script][custom-script-tin] based on
[**RSeqC**](#third-party-software-used).

> TIN is conceptually similar to RIN (RNA integrity number) but provides
> transcript-level measurement of RNA quality and is more sensitive for
> low-quality RNA samples.
>
> - TIN score of a transcript reflects RNA integrity of the transcript
> - Median TIN score across all transcripts reflects RNA integrity of library
> - TIN ranges from 0 (the worst) to 100 (the best)
> - A TIN of 60 means that 60% of the transcript would have been covered if
>   the read coverage were uniform
> - A TIN of 0 will be assigned if the transcript has no coverage or coverage
>   is below threshold

- **Input**
  - Alignments file (`.bam`); from [**map_genome_star**](#map_genome_star)
  - BAM index file (`.bam.bai`); from
    [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
  - Transcript annotations file (12-column `.bed`); from
    [**extract_transcripts_as_bed12**](#extract_transcripts_as_bed12)
- **Output**
  - TIN score table (custom `tsv`); used in
    [**merge_TIN_scores**](#merge_tin_scores)

#### `merge_TIN_scores`

Merges TIN score tables for all samples with [custom script][custom-script-tin].

- **Input**
  - TIN score table (custom `tsv`); per sample; from
    [**calculate_TIN_scores**](#calculate_tin_scores)
- **Output**
  - TIN score table (custom `tsv`); for all samples; used in
    [**plot_TIN_scores**](#plot_tin_scores)

#### `plot_TIN_scores`

Generate sample-wise [box plots](https://en.wikipedia.org/wiki/Box_plot) of
TIN scores with [custom script][custom-script-tin].

- **Input**
  - TIN score table (custom `tsv`); for all samples; from
    [**merge_TIN_scores**](#merge_tin_scores)
- **Output**
  - TIN score box plots (`.pdf` and `.png`); used in
    [**multiqc_report**](#multiqc_report)

#### `salmon_quantmerge_genes`

Merge gene-level expression estimates for all samples with
[**Salmon**](#third-party-software-used).

> Rule is run once per sequencing mode

- **Input**
  - Gene expression tables (custom `.tsv`) for samples of same sequencing mode;
    from [**quantification_salmon**](#quantification_salmon)
- **Output**
  - Gene TPM table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)
  - Gene read count table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)

#### `salmon_quantmerge_transcripts`

Merge transcript-level expression estimates for all samples with
[**Salmon**](#third-party-software-used).

> Rule is run once per sequencing mode

- **Input**
  - Transcript expression tables (custom `.tsv`) for samples of same sequencing
    mode; from [**quantification_salmon**](#quantification_salmon)
- **Output**
  - Transcript TPM table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)
  - Transcript read count table (custom `.tsv`); used in
    [**multiqc_report**](#multiqc_report)

#### `generate_alfa_index`

Create index for [**ALFA**](#third-party-software-used).

- **Input**
  - Gene annotation file (`.gtf`)
  - Chromosome length table `chrNameLength.txt`; from
    [**create_index_star**](#create_index_star)
- **Output**
  - ALFA index; stranded and unstranded; used in
    [**alfa_qc**](#alfa_qc)

#### `alfa_qc`

Annotate alignments with [**ALFA**](#third-party-software-used).

> For details on output plots, see [ALFA documentation][docs-alfa].

- **Input**
  - Coverage files, renamed (`.bg`); from
    [**rename_star_rpm_for_alfa**](#rename_star_rpm_for_alfa)
  - ALFA index, stranded; from [**generate_alfa_index**](#generate_alfa_index)
- **Parameters**
  - `-s`: library orientation; specified by user in sample table column
    `kallisto_directionality`
- **Output**
  - Figures for biotypes and feature categories (`.pdf`)
  - Feature counts table (custom `.tsv`); used in
    [**alfa_qc_all_samples**](#alfa_qc_all_samples)

#### `alfa_qc_all_samples`

Combines output of all samples with [**ALFA**](#third-party-software-used).

- **Input**
  - Feature counts table (custom `.tsv`); from [**alfa_qc**](#alfa_qc)
- **Output**
  - Figures for biotypes and feature categories (`.pdf`); summarized for all
    samples together; used in [**alfa_concat_results**](#alfa_concat_results)

#### `alfa_concat_results`

Concatenate and convert ALFA output plots into single plot with
[**ImageMagick**](#third-party-software-used).

- **Input**
  - Figures for biotypes and feature categories (`.pdf`); for individual and
    summarized for all samples
- **Output**
  - ALFA plot (`.png`), combined; used in [**multiqc_report**](#multiqc_report)

#### `prepare_multiqc_config`

Prepare config file for [**MultiQC**](#third-party-software-used).

> Local rule

- **Input**
  - Directories created during
    [**prepare_files_for_report**](#prepare_files_for_report)
- **Output**
  - Config file (`.yaml`); used in [**multiqc_report**](#multiqc_report)

#### `multiqc_report`

Prepare interactive report from results and logs with
[**MultiQC**](#third-party-software-used).

- **Input**
  - Config file (`.yaml`); from
    [**prepare_multiqc_config**](#prepare_multiqc_config)
  - ALFA plot, combined (`.png`); from
    [**alfa_concat_results**](#alfa_concat_results)
  - Coverage file (`.bg`); from [**star_rpm**](#star_rpm)
  - FastQC output directory with report (`.txt`) and figures (`.png`); from
    [**fastqc**](#fastqc)
  - Gene TPM table (custom `.tsv`); from
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Gene read count table (custom `.tsv`); from
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Pseudoalignments file (`.sam`); from
    [**genome_quantification_kallisto**](#genome_quantification_kallisto)
  - TIN score box plots (`.pdf` and `.png`); from
    [**plot_TIN_scores**](#plot_tin_scores)
  - Transcript TPM table (custom `.tsv`); from
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)
  - Transcript read count table (custom `.tsv`); from
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)

- **Output**
  - Directory with automatically generated `.html` report

#### `finish`

Target rule as required by [Snakemake][docs-snakemake-target-rule].

> Local rule

- **Input**
  - MultiQC report; from [**multiqc_report**](#multiqc_report)
  - Coverage files, one per strand and sample (`.bw`); used in
    [**prepare_bigWig**](#prepare_bigwig)

### Sequencing mode-specific

> Steps described here have two variants, one with the specified names for
> samples prepared with a single-end sequencing protocol, one with `pe_`
> prepended to the specified name for samples prepared with a paired-end
> sequencing protocol.

#### `remove_adapters_cutadapt`

Remove adapter sequences from reads with
[**Cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from [**start**](#start)
- **Parameters**
  - Adapters to be removed; specify in sample table columns `fq1_3p`, `fq1_5p`,
    `fq2_3p`, `fq2_5p`
- **Output**
  - Reads file (`.fastq.gz`); used in
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
- **Non-configurable & non-default**
  - `-e 0.1`: maximum error-rate of 10%
  - `-j 8`: use 8 threads
  - `-m 10`: Discard processed reads that are shorter than 10
  - `-n 2`: search for all the given adapter sequences repeatedly, either until
    no adapter match was found or until 2 rounds have been performed.
  - `--pair-filter=any`: **(paired-end only)** filtering criteria must apply to
    any of the two reads in order for a read pair to be discarded

#### `remove_polya_cutadapt`

Remove poly(A) tails from reads with
[**Cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_adapters_cutadapt**](#remove_adapters_cutadapt)
- **Parameters**
  - Poly(A) stretches to be removed; specify in sample table columns `fq1_polya`
    and `fq2_polya`
- **Output**
  - Reads file (`.fastq.gz`); used in
    [**genome_quantification_kallisto**](#genome_quantification_kallisto),
    [**map_genome_star**](#map_genome_star) and
    [**quantification_salmon**](#quantification_salmon)
- **Non-configurable & non-default**
  - `-e 0.1`: maximum error-rate of 10%
  - `-j 8`: use 8 threads
  - `-m 10`: Discard processed reads that are shorter than 10
  - `-n 1`: search for all the given adapter sequences repeatedly, either until
    no adapter match was found or until 1 round has been performed.
  - `--pair-filter=any`: **(paired-end only)** filtering criteria must apply to
    any of the two reads in order for a read pair to be discarded
  - `-O 1`: **(single-end only)** minimal overlap of 1

#### `map_genome_star`

Align short reads to reference genome and/or transcriptome with
[**STAR**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Index; from [**create_index_star**](#create_index_star)
- **Parameters**
  - `--outFilterMultimapNmax`: maximum number of multiple alignments allowed;
    if exceeded, read is considered unmapped; specify in sample table column
    `multimappers`
  - `--alignEndsType`: one of `Local` (standard local alignment with
    soft-clipping allowed) or `EndToEnd` (force end-to-end read alignment, do
    not soft-clip); specify in sample table column `soft_clip`
  - `--twopassMode`: one of `None` (1-pass mapping) or `Basic` (basic 2-pass
    mapping, with all 1st-pass junctions inserted into the genome indices on
    the fly); specify in sample table column `pass_mode`
- **Output**
  - Aligned reads file (`.bam`); used in
    [**calculate_TIN_scores**](#calculate_TIN_scores),
    [**index_genomic_alignment_samtools**](#index_genomic_alignment_samtools)
    and [**star_rpm**](#star_rpm)
  - STAR log file
- **Non-configurable & non-default**
  - `--outSAMunmapped=None`: do not output unmapped reads in SAM file
  - `--outFilterMultimapScoreRange=0`: the score range below the maximum score
    for multimapping alignments
  - `--outSAMattributes=All`: NH HI AS nM NM MD jM jI MC ch
  - `--outStd=BAM_SortedByCoordinate`: which output will be directed to `STDOUT`
  - `--outSAMtype=BAM_SortedByCoordinate`: type of SAM/BAM output
  - `--outFilterMismatchNoverLmax=0.04`: alignment will be output only if its
    ratio of mismatches to *mapped* length is less than or equal to this value
  - `--outFilterScoreMinOverLread=0.3`: same as outFilterScoreMin, but
    normalized to read length (sum of mates’ lengths for paired-end reads)
  - `--outFilterMatchNminOverLread=0.3`: minimal fraction of aligned bases
  - `--outFilterType=BySJout`: reduces the number of ”spurious” junctions
  - `--outReadsUnmapped=None`: do not output unmapped reads

#### `quantification_salmon`

Estimate transcript- and gene-level expression with
[**Salmon**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Filtered annotation file (`.gtf`)
  - Index; from [**create_index_salmon**](#create_index_salmon)
- **Parameters**
  - `libType`: see [Salmon manual][docs-salmon] for allowed values; specify in
    sample table column `libtype`
  - `--fldMean`: mean of distribution of fragment lengths; specify in sample
    table column `mean` **(single-end only)**
  - `--fldSD`: standard deviation of distribution of fragment lengths; specify
    in sample table column `sd` **(single-end only)**
- **Output**
  - Gene expression table (custom `.tsv`); used in
    [**salmon_quantmerge_genes**](#salmon_quantmerge_genes)
  - Transcript expression table (custom `.tsv`); used in
    [**salmon_quantmerge_transcripts**](#salmon_quantmerge_transcripts)
- **Non-configurable & non-default**
  - `--seqBias`: [correct for sequence specific
    biases](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias)
  - `--validateMappings`: enables selective alignment of the sequencing reads
    when mapping them to the transcriptome; this can improve both the
    sensitivity and specificity of mapping and, as a result, can [improve
    quantification
    accuracy](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings).
  - `--writeUnmappedNames`: write out the names of reads (or mates in paired-end
    reads) that do not map to the transcriptome.

#### `genome_quantification_kallisto`

Generate pseudoalignments of reads to transcripts with
[**kallisto**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq.gz`); from
    [**remove_polya_cutadapt**](#remove_polya_cutadapt)
  - Index; from [**create_index_kallisto**](#create_index_kallisto)
- **Parameters**
  - `directionality`; specify in sample table column `kallisto_directionality`
  - `-l`: mean of distribution of fragment lengths; specify in sample table
    column `mean` **(single-end only)**
  - `-s`: standard deviation of distribution of fragment lengths; specify in
    sample table column `sd` **(single-end only)**
- **Output**
  - Pseudoalignments file (`.sam`); used in
    [**multiqc_report**](#multiqc_report)

[code-alfa]: <https://github.com/biocompibens/ALFA>
[code-bedgraphtobigwig]: <https://github.com/ucscGenomeBrowser/kent>
[code-bedtools]: <https://github.com/arq5x/bedtools2>
[code-cutadapt]: <https://github.com/marcelm/cutadapt>
[code-gffread]: <https://github.com/gpertea/gffread>
[code-fastqc]: <https://github.com/s-andrews/FastQC>
[code-imagemagick]: <https://github.com/ImageMagick/ImageMagick/>
[code-kallisto]: <https://github.com/pachterlab/kallisto>
[code-multiqc]: <https://github.com/ewels/MultiQC>
[code-rseqc]: <http://rseqc.sourceforge.net/>
[code-salmon]: <https://github.com/COMBINE-lab/salmon>
[code-samtools]: <https://github.com/samtools/samtools>
[code-star]: <https://github.com/alexdobin/STAR>
[custom-script-gtf-to-bed12]: <https://git.scicore.unibas.ch/zavolan_group/tools/gtf_transcript_type_to_bed12>
[custom-script-tin]: <https://git.scicore.unibas.ch/zavolan_group/tools/tin_score_calculation>
[docs-alfa]: <https://github.com/biocompibens/ALFA#manual>
[docs-bedgraphtobigwig]: <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>
[docs-bedtools]: <https://bedtools.readthedocs.io/en/latest/>
[docs-cutadapt]: <https://cutadapt.readthedocs.io/en/stable/>
[docs-gffread]: <http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread>
[docs-fastqc]: <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/>
[docs-imagemagick]: <https://imagemagick.org/>
[docs-kallisto]: <http://pachterlab.github.io/kallisto/manual.html>
[docs-multiqc]: <https://multiqc.info/docs/>
[docs-rseqc]: <http://rseqc.sourceforge.net/#usage-information>
[docs-salmon]: <https://salmon.readthedocs.io/en/latest/>
[docs-salmon-selective-alignment]: <https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/>
[docs-samtools]: <http://www.htslib.org/doc/samtools.html>
[docs-snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[docs-snakemake-target-rule]: <https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-7-adding-a-target-rule>
[docs-star]: <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
[docs-star-rpm-norm]: <https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html>
[license-bsd2]: <https://opensource.org/licenses/BSD-2-Clause>
[license-gpl2]: <https://opensource.org/licenses/GPL-2.0>
[license-gpl3]: <https://opensource.org/licenses/GPL-3.0>
[license-imagemagick]: <https://github.com/ImageMagick/ImageMagick/blob/master/LICENSE>
[license-mit]: <https://opensource.org/licenses/MIT>
[pub-alfa]: <https://doi.org/10.1186/s12864-019-5624-2>
[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-kallisto]: <https://doi.org/10.1038/nbt.3519>
[pub-multiqc]: <https://doi.org/10.1093/bioinformatics/btw354>
[pub-rseqc]: <https://doi.org/10.1093/bioinformatics/bts356>
[pub-salmon]: <https://doi.org/10.1038/nmeth.4197>
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-star]: <https://doi.org/10.1093/bioinformatics/bts635>
[rule-graph]: images/rule_graph.svg
