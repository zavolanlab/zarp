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
dependency for ZARP is the [Mamba][mamba] package manager (version 1), which needs to be installed in
the `base` conda environment with:

```bash
conda install mamba=1 -n base -c conda-forge
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

## Running ZARP without ZARP-cli

You can also trigger ZARP without ZARP-cli. This is convenient for users who have some experience with snakemake and don't want to use a CLI to trigger their runs. Please head over to the [ZARP](https://zavolanlab.github.io/zarp/) documentation to learn how to start ZARP.

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
