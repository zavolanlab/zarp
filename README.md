[![ci](https://github.com/zavolanlab/zarp/workflows/CI/badge.svg?branch=dev)](https://github.com/zavolanlab/zarp/actions?query=workflow%3Aci)
[![GitHub license](https://img.shields.io/github/license/zavolanlab/zarp?color=orange)](https://github.com/zavolanlab/zarp/blob/dev/LICENSE)
[![DOI:biorxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2021.11.18.469017-informational)](https://doi.org/10.1101/2021.11.18.469017)
[![DOI:zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.5703358-informational)](https://doi.org/10.5281/zenodo.5703358)
[![DOI:workflowhub](https://img.shields.io/badge/WorkflowHub-10.48546%2Fworkflowhub.workflow.447.1-informational)](https://doi.org/10.48546/workflowhub.workflow.447.1)

<div align="left">
    <img width="20%" align="left" src=images/zarp_logo.svg>
</div> 

**ZARP** ([Zavolab][zavolan-lab] Automated RNA-seq Pipeline) is a generic
RNA-Seq analysis workflow that allows users to process and analyze Illumina
short-read sequencing libraries with minimum effort. Better yet: With our
companion [**ZARP-cli**](https://github.com/zavolanlab/zarp-cli) command line
interface, you can start ZARP runs with the simplest and most intuitive
commands.

_RNA-seq analysis doesn't get simpler than that!_

ZARP relies on publicly available bioinformatics tools and currently handles
single or paired-end stranded bulk RNA-seq data. The workflow is developed in
[Snakemake][snakemake], a widely used workflow management system in the
bioinformatics community.

ZARP will pre-process, align and quantify your single- or paired-end stranded
bulk RNA-seq sequencing libraries with publicly available state-of-the-art
bioinformatics tools. ZARP's browser-based rich reports and visualitations will
give you meaningful initial insights in the quality and composition of your
sequencing experiments - fast and simple. Whether you are an experimentalist
struggling with large scale data analysis or an experienced bioinformatician,
when there's RNA-seq data to analyze, just _zarp 'em_!

<div align="center">
    <img width="60%" src=images/zarp_schema.png>
</div> 

> **Note:** For a more detailed description of each step, please refer to the [workflow
> documentation][pipeline-documentation].

# Requirements

The workflow has been tested on:
- CentOS 7.5
- Debian 10
- Ubuntu 16.04, 18.04

> **NOTE:**
> Currently, we only support **Linux** execution. 


# Installation

> **IMPORTANT: Rather than installing the ZARP workflow as described in this section, we
> recommend installing [ZARP-cli](https://github.com/zavolanlab/zarp-cli) for most use
> cases!** If you follow its [installation
> instructions](https://zavolanlab.github.io/zarp-cli/guides/installation/), you can
> skip the instructions below.

## 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone https://github.com/zavolanlab/zarp.git
cd zarp
```

## 2. Conda and Mamba installation

Workflow dependencies can be conveniently installed with the [Conda][conda]
package manager. We recommend that you install [Miniconda][miniconda-installation] 
for your system (Linux). Be sure to select Python 3 option. 
The workflow was built and tested with `miniconda 4.7.12`.
Other versions are not guaranteed to work as expected.

Given that Miniconda has been installed and is available in the current shell the first
dependency for ZARP is the [Mamba][mamba] package manager, which needs to be installed in
the `base` conda environment with:

```bash
conda install mamba -n base -c conda-forge
```

## 3. Dependencies installation

For improved reproducibility and reusability of the workflow,
each individual step of the workflow runs either in its own [Singularity][singularity]
container or in its own [Conda][conda] virtual environemnt. 
As a consequence, running this workflow has very few individual dependencies. 
The **container execution** requires Singularity to be installed on the system where the workflow is executed. 
As the functional installation of Singularity requires root privileges, and Conda currently only provides Singularity
for Linux architectures, the installation instructions are slightly different depending on your system/setup:

### For most users

If you do *not* have root privileges on the machine you want
to run the workflow on *or* if you do not have a Linux machine, please [install
Singularity][singularity-install] separately and in privileged mode, depending
on your system. You may have to ask an authorized person (e.g., a systems
administrator) to do that. This will almost certainly be required if you want
to run the workflow on a high-performance computing (HPC) cluster. 

> **NOTE:**
> The workflow has been tested with the following Singularity versions:  
>  * `v2.6.2`
>  * `v3.5.2`

After installing Singularity, install the remaining dependencies with:
```bash
mamba env create -f install/environment.yml
```


### As root user on Linux

If you have a Linux machine, as well as root privileges, (e.g., if you plan to
run the workflow on your own computer), you can execute the following command
to include Singularity in the Conda environment:

```bash
mamba env update -f install/environment.root.yml
```

## 4. Activate environment

Activate the Conda environment with:

```bash
conda activate zarp
```

# Extra installation steps (optional)

## 5. Non-essential dependencies installation

Most tests have additional dependencies. If you are planning to run tests, you
will need to install these by executing the following command _in your active
Conda environment_:

```bash
mamba env update -f install/environment.dev.yml
```

## 6. Successful installation tests

We have prepared several tests to check the integrity of the workflow and its
components. These can be found in subdirectories of the `tests/` directory. 
The most critical of these tests enable you to execute the entire workflow on a 
set of small example input files. Note that for this and other tests to complete
successfully, [additional dependencies](#installing-non-essential-dependencies) 
need to be installed. 
Execute one of the following commands to run the test workflow 
on your local machine:
* Test workflow on local machine with **Singularity**:
```bash
bash tests/test_integration_workflow/test.local.sh
```
* Test workflow on local machine with **Conda**:
```bash
bash tests/test_integration_workflow_with_conda/test.local.sh
```
Execute one of the following commands to run the test workflow 
on a [Slurm][slurm]-managed high-performance computing (HPC) cluster:

* Test workflow with **Singularity**:

```bash
bash tests/test_integration_workflow/test.slurm.sh
```
* Test workflow with **Conda**:

```bash
bash tests/test_integration_workflow_with_conda/test.slurm.sh
```

> **NOTE:** Depending on the configuration of your Slurm installation you may
> need to adapt file `slurm-config.json` (located directly under `profiles`
> directory) and the arguments to options `--cores` and `--jobs`
> in the file `config.yaml` of a respective profile.
> Consult the manual of your workload manager as well as the section of the
> Snakemake manual dealing with [profiles].

# Running the workflow on your own samples

## Running ZARP with ZARP-cli (recommended)

Head over to the [ZARP-cli](https://zavolanlab.github.io/zarp-cli/) to learn how to
start ZARP runs with very simple commands, like:

```sh
zarp SRR23590181
```

## Running ZARP without ZARP-cli

1. Assuming that your current directory is the workflow repository's root directory,
create a directory for your workflow run and move into it with:

    ```bash
    mkdir config/my_run
    cd config/my_run
    ```

2. Create an empty sample table and a workflow configuration file:

    ```bash
    touch samples.tsv
    touch config.yaml
    ```

3. Use your editor of choice to populate these files with appropriate
values. Have a look at the examples in the `tests/` directory to see what the
files should look like, specifically:

    - [samples.tsv](tests/input_files/samples.tsv)
    - [config.yaml](tests/input_files/config.yaml)

    - For more details and explanations, refer to the [pipeline-documentation]


4. Create a runner script. Pick one of the following choices for either local
or cluster execution. Before execution of the respective command, you need to
remember to update the argument of the `--singularity-args` option of a
respective profile (file: `profiles/{profile}/config.yaml`) so that
it contains a comma-separated list of _all_ directories
containing input data files (samples and any annotation files etc) required for
your run.

    Runner script for _local execution_:

    ```bash
    cat << "EOF" > run.sh
    #!/bin/bash

    snakemake \
        --profile="../../profiles/local-singularity" \
        --configfile="config.yaml"

    EOF
    ```

    **OR**

    Runner script for _Slurm cluster exection_ (note that you may need
    to modify the arguments to `--jobs` and `--cores` in the file:
    `profiles/slurm-singularity/config.yaml` depending on your HPC
    and workload manager configuration):

    ```bash
    cat << "EOF" > run.sh
    #!/bin/bash
    mkdir -p logs/cluster_log
    snakemake \
        --profile="../profiles/slurm-singularity" \
        --configfile="config.yaml"
    EOF
    ```

    > Note: When running the pipeline with *conda* you should use `local-conda` and
    `slurm-conda` profiles instead.

    > Note: The slurm profiles are adapted to a cluster that uses the quality-of-service (QOS) keyword. If QOS is not supported by your slurm instance, you have to remove all the lines with "qos" in `profiles/slurm-config.json`.

5. Start your workflow run:

    ```bash
    bash run.sh
    ```

## Output files

 After running the ZARP workflow, you will find several output files in the specified output directory. The output directory is defined in the `config.yaml` file and it is normally called `results`. Here are some of the key output files:

- **Alignment files**: These files contain the aligned reads for each sample in BAM format. They provide information about the mapping of the reads to the reference genome.

- **Quantification files**: These files contain the gene and transcript level expression values for each sample. They provide information about the abundance of each gene / transcript in the RNA-seq data.

- **Quality control reports**: ZARP generates comprehensive quality control reports that provide insights into the quality and composition of the sequencing experiments. These reports include metrics such as read quality, alignment statistics, and gene expression summaries.

- **Visualization files**: ZARP also generates visualizations to help you interpret the results. These visualizations include plots, statistics and interactive browser-based reports.

After a run you will find the following structure within the `results` directory:

```bash
.
├── multiqc_config.yaml
└── mus_musculus
    ├── multiqc_summary
    ├── samples
    ├── summary_kallisto
    ├── summary_salmon
    └── zpca
```

A descrpition of the different directories is shown below:

- `results`: The main output directory for the ZARP workflow.
    - `mus_musculus`: A subdirectory for the organism-specific results.
        - `multiqc_summary`: Summary files generated by MultiQC.
        - `samples`: Sample specific outputs. A directory is created for each sample.
        - `summary_kallisto`: Summary files for Kallisto quantifications.
        - `summary_salmon`: Summary files for Salmon quantifications.
        - `zpca`: Output files for ZARP's principal component analysis.

### QC outputs

Within the `multiqc_summary` directory, you will find an interactive HTML file (`multiqc_report.html`) with various QC metrics that can help you interpret your results. An example file is shown below

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc.png>
</div>

On the left you can find a navigation bar that takes you into different sections and subsections of the tools.

- The `General Statistics` section contains a summary of most tools and you can find statistics on mapped reads, percent of duplicate reads, percent of adapters trimmed for various tools.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_general_statistics.png>
</div>

- The `FastQC: raw reads` section contains plots and quality statistics of the fastq files. Some examples are shown below like the number of duplicate reads in an experiment, the average quality of the fastq files per position, or the percent of GC content.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_fastqc_sequence_counts_plot.png>
</div>

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_fastqc_per_base_sequence_quality_plot.png>
</div>

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_fastqc_per_sequence_gc_content_plot.png>
</div>

- The `Cutadapt: adapter removal` and `Cutadapt: polyA tails removal` shows the number or the percentage of the reads trimmed

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_cutadapt_filtered_reads_plot.png>
</div>


- The `FastQC: trimmed reads` section contains plots and quality statistics of the fastq files after adapter trimming. The plots are similar to the section `FastQC: raw reads`.

- The `STAR` section shows the number and percentage of reads that are mapped using the STAR aligner.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_star_alignment_plot.png>
</div>

- The `ALFA` section shows the number of reads mapped to genomic categories (stop codon, 5'-UTR, CDS, intergenic, etc.) and gene biotypes (protein coding genes, miRNA , tRNA, etc.) for unique reads and multimappers.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_alfa_categories.png>
</div>

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_alfa_biotypes.png>
</div>

- The `TIN` section shows the Transcript Integrity Number of the samples.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_tin_score.png>
</div>

- The `Salmon` section shows the fragment length distribution of the reads

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_salmon_plot.png>
</div>

- The `Kallisto` section shows the number of reads that were aligned

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_kallisto_alignment.png>
</div>

- Finally the `zpca` salmon and kallisto sections show PCA plots for expression levels of genes and transcripts.

<div align="center">
    <img width="80%" src=images/output_files/zarp_multiqc_zpca.png>
</div>

### Gene and transcript estimate outputs

Within the `summary_kallisto` directory, you can find the following files:
- `genes_counts.tsv`: Matrix with the gene counts. The first column (index) contains the gene names and the first row (column) contains the sample names. This file can later be used for downstream differential expression analysis. 
- `genes_tpm.tsv`: Matrix with the gene TPM estimates.
- `transcripts_counts.tsv`: Matrix with the transcript counts. The first column (index) contains the transcript names and the first row (column) contains the sample names. This file can later be used for downstream differential transcript analysis.
- `transcripts_tpm.tsv`: Matrix with the transcript TPM estimates.
- `tx2geneID.tsv`: A table mapping transcript IDs to gene IDs.

Within the `summary_salmon/quantmerge` directory, you can find the following files:
- `genes_numreads.tsv`: Matrix with the gene counts. The first column (index) contains the gene names and the first row (column) contains the sample names. This file can later be used for downstream differential expression analysis. 
- `genes_tpm.tsv`: Matrix with the gene TPM estimates. 
- `transcripts_numreads.tsv`: Matrix with the transcript counts. The first column (index) contains the transcript names and the first row (column) contains the sample names. This file can later be used for downstream differential transcript analysis.
- `transcripts_tpm.tsv`: Matrix with the transcript TPM estimates.

### Alignment outputs

Within the `samples` directory, you can find a directory for each sample, and within these directories you can find the output files of the individual steps. Some alignment files can be easily used to open in a genome browser for other downstream analysis:
- In the `map_genome` directory you can find a file with the suffix `.Aligned.sortedByCoord.out.bam` and the corresponding indexed (`.bai`) file. This is the output of the STAR aligner. 
- In the `bigWig` directory you can find two folders. `UniqueMappers` and `MultimappersIncluded`. Within these files you find the bigWig files for the plus and minus strand. These files are convenient to load in a genome browser (like igv) to view the genome coverage of the mappings.


# Sample downloads from SRA

An independent Snakemake workflow `workflow/rules/sra_download.smk` is included
for the download of sequencing libraries from the Sequence Read Archive and
conversion into FASTQ.

The workflow expects the following parameters in the configuration file:
* `samples`, a sample table (tsv) with column *sample* containing *SRR*
  identifiers (ERR and DRR are also supported), see
  [example](tests/input_files/sra_samples.tsv).
* `outdir`, an output directory
* `samples_out`, a pointer to a modified sample table with the locations of
  the corresponding FASTQ files
* `cluster_log_dir`, the cluster log directory.

For executing the example with Conda environments, one can use the following
command (from within the activated `zarp` Conda environment):

```bash
snakemake --snakefile="workflow/rules/sra_download.smk" \
          --profile="profiles/local-conda" \
          --config samples="tests/input_files/sra_samples.tsv" \
                   outdir="results/sra_downloads" \
                   samples_out="results/sra_downloads/sra_samples.out.tsv" \
                   log_dir="logs" \
                   cluster_log_dir="logs/cluster_log"
```

Alternatively, change the argument to `--profile` from `local-conda` to
`local-singularity` to execute the workflow steps within Singularity
containers.

After successful execution, `results/sra_downloads/sra_samples.out.tsv` should
contain:

```tsv
sample  fq1     fq2
SRR18552868     results/sra_downloads/compress/SRR18552868/SRR18552868.fastq.gz 
SRR18549672     results/sra_downloads/compress/SRR18549672/SRR18549672_1.fastq.gz       results/sra_downloads/compress/SRR18549672/SRR18549672_2.fastq.gz
ERR2248142      results/sra_downloads/compress/ERR2248142/ERR2248142.fastq.gz 
```


# Metadata completion with HTSinfer
An independent Snakemake workflow `workflow/rules/htsinfer.smk` that populates the `samples.tsv` required by ZARP with the sample specific parameters `seqmode`, `f1_3p`, `f2_3p`, `organism`, `libtype` and `index_size`. Those parameters are inferred from the provided `fastq.gz` files by [HTSinfer][hts-infer].

> Note: The workflow uses the implicit temporary directory 
from snakemake, which is called with [resources.tmpdir].


The workflow expects the following config:
* `samples`, a sample table (tsv) with column *sample* containing sample identifiers, as well as columns *fq1* and *fq2* containing the paths to the input fastq files
see example [here](tests/input_files/sra_samples.tsv). If the table contains further ZARP compatible columns (see [pipeline documentation][sample-doc]), the values specified there by the user are given priority over htsinfer's results. 
* `outdir`, an output directory
* `samples_out`, path to a modified sample table with inferred parameters
* `records`, set to 100000 per default
  
For executing the example one can use the following
(with activated *zarp* environment):
```bash
cd tests/test_htsinfer_workflow
snakemake \
    --snakefile="../../workflow/rules/htsinfer.smk" \
    --restart-times=0 \
    --profile="../../profiles/local-singularity" \
    --config outdir="results" \
             samples="../input_files/htsinfer_samples.tsv" \
             samples_out="samples_htsinfer.tsv" \
    --notemp \
    --keep-incomplete
```

However, this call will exit with an error, as not all parameters can be inferred from the example files. The argument `--keep-incomplete` makes sure the `samples_htsinfer.tsv` file can nevertheless be inspected. 

After successful execution - if all parameters could be either inferred or were specified by the user - `[OUTDIR]/[SAMPLES_OUT]` should contain a populated table with parameters `seqmode`, `f1_3p`, `f2_3p`, `organism`, `libtype` and `index_size` for all input samples as described in the [pipeline documentation][sample-doc].



[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[hts-infer]: <https://github.com/zavolanlab/htsinfer>
[profiles]: <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>
[mamba]: <https://github.com/mamba-org/mamba>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[rule-graph]: images/rule_graph.svg
[zarp-logo]: images/zarp_logo.svg
[zarp-schema]: images/zarp_schema.svg
[sample-doc]: pipeline_documentation.md#read-sample-table
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[singularity]: <https://sylabs.io/singularity/>
[singularity-install]: <https://sylabs.io/guides/3.5/admin-guide/installation.html>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[zavolan-lab]: <https://www.biozentrum.unibas.ch/research/researchgroups/overview/unit/zavolan/research-group-mihaela-zavolan/>
[pipeline-documentation]: pipeline_documentation.md
[resources.tmpdir]: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?#standard-resources>
