[![ci](https://github.com/zavolanlab/zarp/workflows/CI/badge.svg?branch=dev)](https://github.com/zavolanlab/zarp/actions?query=workflow%3Aci)
[![GitHub license](https://img.shields.io/github/license/zavolanlab/zarp?color=orange)](https://github.com/zavolanlab/zarp/blob/dev/LICENSE)
[![Static Badge](https://img.shields.io/badge/f1000-10.12688/f1000research.149237.1-blue)](https://doi.org/10.12688/f1000research.149237.1)
[![DOI:zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.10797025-informational)](https://doi.org/10.5281/zenodo.10797025)
[![DOI:workflowhub](https://img.shields.io/badge/WorkflowHub-10.48546%2Fworkflowhub.workflow.447.1-informational)](https://doi.org/10.48546/workflowhub.workflow.447.1)

<div align="left">
    <img width="20%" align="left" src=images/zarp_logo.svg>
</div>

**ZARP** ([Zavolab][zavolan-lab] Automated RNA-seq Pipeline) is a generic RNA-Seq analysis workflow that allows users to process and analyze Illumina short-read sequencing libraries with minimum effort. Better yet: With our companion [**ZARP-cli**](https://github.com/zavolanlab/zarp-cli) command line interface, you can start ZARP runs with the simplest and most intuitive commands.

_RNA-seq analysis doesn't get simpler than that!_

ZARP relies on publicly available bioinformatics tools and currently handles single or paired-end stranded bulk RNA-seq data. The workflow is developed in [Snakemake][snakemake], a widely used workflow management system in the bioinformatics community.

ZARP will pre-process, align and quantify your single- or paired-end stranded bulk RNA-seq sequencing libraries with publicly available state-of-the-art bioinformatics tools. ZARP's browser-based rich reports and visualitations will give you meaningful initial insights in the quality and composition of your sequencing experiments - fast and simple. Whether you are an experimentalist struggling with large scale data analysis or an experienced bioinformatician, when there's RNA-seq data to analyze, just _ZARP 'em_!

<div align="center">
    <img width="60%" src=images/zarp_schema.png>
</div> 

# Documentation

For the full documentation please visit the [ZARP website](https://zavolanlab.github.io/zarp).

# Quick installation

> **IMPORTANT: Rather than installing the ZARP workflow as described in this section, we
> recommend installing [ZARP-cli](https://github.com/zavolanlab/zarp-cli) for most use
> cases!** If you follow its [installation
> instructions](https://zavolanlab.github.io/zarp-cli/guides/installation/), you can
> skip the instructions below.

Quick installation requires the following:
- Linux
- Git
- [Conda][conda] >= 22.11.1
- [Mamba][mamba] >=1.3.0 <2
- [Singularity][singularity] >=3.5.2  (Required only if you want to use Singulaarity for the dependencies)

```bash
git clone https://github.com/zavolanlab/zarp.git
cd zarp
mamba env create -f install/environment.yml
conda activate zarp
```

# Basic usage

You can trigger ZARP without ZARP-cli. This is convenient for users who have some experience with Snakemake and don't want to use a CLI to trigger their runs. Extensive documentation of the usage is available in the [usage documentation](https://zavolanlab.github.io/zarp/guides/usage/), while below you can find the basic steps to trigger a run.

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

    - [samples.tsv](https://github.com/zavolanlab/zarp/blob/dev/tests/input_files/samples.tsv)
    - [config.yaml](https://github.com/zavolanlab/zarp/blob/dev/tests/input_files/config.yaml)


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

    > Note: When running the pipeline with *Conda* you should use `local-conda` and
    `slurm-conda` profiles instead.

    > Note: The slurm profiles are adapted to a cluster that uses the quality-of-service (QOS) keyword. If QOS is not supported by your slurm instance, you have to remove all the lines with "qos" in `profiles/slurm-config.json`.

5. Start your workflow run:

    ```bash
    bash run.sh
    ```

## Contributing

This project lives off your contributions, be it in the form of bug reports,
feature requests, discussions, or fixes and other code changes. Please refer
to the [contributing guidelines](CONTRIBUTING.md) if you are interested to
contribute. Please mind the [code of conduct](CODE_OF_CONDUCT.md) for all
interactions with the community.

## Contact

For questions or suggestions regarding the code, please use the
[issue tracker][issue-tracker]. For any other inquiries, please contact us
by [email][contact].

&copy; 2021 [Zavolab, Biozentrum, University of Basel][zavolab]


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
[zavolab]: <https://www.biozentrum.unibas.ch/research/researchgroups/overview/unit/zavolan/research-group-mihaela-zavolan/>
[contact]: <mailto:zavolab-biozentrum@unibas.ch>