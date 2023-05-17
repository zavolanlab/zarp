[![ci](https://github.com/zavolanlab/zarp/workflows/CI/badge.svg?branch=dev)](https://github.com/zavolanlab/zarp/actions?query=workflow%3Aci)
[![GitHub license](https://img.shields.io/github/license/zavolanlab/zarp?color=orange)](https://github.com/zavolanlab/zarp/blob/dev/LICENSE)
[![DOI:biorxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2021.11.18.469017-informational)](https://doi.org/10.1101/2021.11.18.469017)
[![DOI:zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.5703358-informational)](https://doi.org/10.5281/zenodo.5703358)
[![DOI:workflowhub](https://img.shields.io/badge/WorkflowHub-10.48546%2Fworkflowhub.workflow.447.1-informational)](https://doi.org/10.48546/workflowhub.workflow.447.1)

<div align="left">
    <img width="20%" align="left" src=images/zarp_logo.svg>
</div> 

**ZARP** ([Zavolan-Lab][zavolan-lab] Automated RNA-Seq Pipeline) is a generic RNA-Seq analysis workflow that allows 
users to process and analyze Illumina short-read sequencing libraries with minimum effort. The workflow relies on 
publicly available bioinformatics tools and currently handles single or paired-end stranded bulk RNA-seq data.
The workflow is developed in [Snakemake][snakemake], a widely used workflow management system in the bioinformatics
community.

According to the current ZARP implementation, reads are analyzed (pre-processed, aligned, quantified) with state-of-the-art
tools to give meaningful initial insights into the quality and composition of an RNA-Seq library, reducing hands-on time for bioinformaticians and giving experimentalists the possibility to rapidly assess their data. Additional reports summarise the results of the individual steps and provide useful visualisations.

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

### Command-line interface dependencies
> A command-line interface, [ZARP-cli](https://github.com/zavolanlab/zarp-cli) is now availble and all you have to do to make use of it is update your environment:

```bash
mamba env update -f install/environment.cli.yml
```

More info at: [NEW: Running ZARP gets even easier!](#new-running-zarp-gets-even-easier)



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

## NEW: Running ZARP gets even easier!
>  You can now run your ZARP workflow using its' new command-line interface, [ZARP-cli](https://github.com/zavolanlab/zarp-cli). To enable ZARP-cli, simply install the additional dependencies by executing the following command _in your active Conda environment_:

```bash
mamba env update -f install/environment.cli.yml
```
Run:
```
zarp --help
```
to see the useful options available, or take a look at [ZARP-cli documentation](https://github.com/zavolanlab/zarp-cli)


## Running ZARP without ZARP-cli

1. Assuming that your current directory is the repository's root directory,
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

# Sample downloads from SRA

An independent Snakemake workflow `workflow/rules/sra_download.smk` is included
for the download of SRA samples with [sra-tools].

> Note: as of Snakemake 7.3.1, only profile conda is supported. 
> Singularity fails because the *sra-tools* Docker container only has `sh` 
but `bash` is required.

> Note: The workflow uses the implicit temporary directory 
from snakemake, which is called with [resources.tmpdir].

The workflow expects the following config:
* `samples`, a sample table (tsv) with column *sample* containing *SRR* identifiers,
see example [here](tests/input_files/sra_samples.tsv).
* `outdir`, an output directory
* `samples_out`, a pointer to a modified sample table with location of fastq files
* `cluster_log_dir`, the cluster log directory.

For executing the example one can use the following
(with activated *zarp* environment):

```bash
snakemake --snakefile="workflow/rules/sra_download.smk" \
          --profile="profiles/local-conda" \
          --config samples="tests/input_files/sra_samples.tsv" \
                   outdir="results/sra_downloads" \
                   samples_out="results/sra_downloads/sra_samples.out.tsv" \
                   log_dir="logs" \
                   cluster_log_dir="logs/cluster_log"
```
After successful execution, `results/sra_downloads/sra_samples.out.tsv` should contain:
```tsv
sample	fq1	fq2
SRR18552868	results/sra_downloads/SRR18552868/SRR18552868.fastq.gz	
SRR18549672	results/sra_downloads/SRR18549672/SRR18549672_1.fastq.gz	results/sra_downloads/SRR18549672/SRR18549672_2.fastq.gz
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
[sra-tools]: <https://github.com/ncbi/sra-tools>
[resources.tmpdir]: <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?#standard-resources>
