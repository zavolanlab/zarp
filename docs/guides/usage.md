# Execution of workflows

This section describes how to run _ZARP_, as well as the two auxiliary
workflow for fetching samples from the Sequence Read Archive and populating a
sparse sample table with inferred sample metadata that are also packaged
within this repository for your convenience.

Finally, we describe how you can use Docker to try to run _ZARP_ on systems
that are not natively supported (e.g., Mac OS).

!!! info "Prerequisites"

    All usage examples in this section assume that you have already
    [installed](./installation.md) _ZARP_.

## How to analyze your RNA-Seq samples?

1. Assuming that your current directory is the workflow repository's root
   directory, create a directory for your workflow run and traverse into it
   with:

    ```bash
    mkdir config/my_run/
    cd config/my_run/
    ```

2. Create an empty sample table and a workflow configuration file:

    ```bash
    touch samples.tsv
    touch config.yaml
    ```

3. Use your editor of choice to populate these files with appropriate
   values. Have a look at the examples in the `tests/` directory to see what
   the files should look like, specifically:

    - [`samples.tsv`][sample-table]
    - [`config.yaml`][config-file]

4. Create a runner script. Pick one of the following choices for either local
   or cluster execution. Before execution of the respective command, you need
   to remember to update the argument of the `--apptainer-args` option of a
   respective profile (file: `profiles/{profile}/config.yaml`) so that
   it contains a comma-separated list of _all_ directories containing input
   data files (samples and any annotation files etc) required for your run.

    Runner script for _local execution_:

    ```bash
    cat << "EOF" > run.sh
    #!/bin/bash

    snakemake \
        --profile="../../profiles/local-apptainer" \
        --configfile="config.yaml"

    EOF
    ```

    **OR**

    Runner script for _Slurm cluster execution_ (note that you may need to
    modify the arguments to `--jobs` and `--cores` in the file
    `profiles/slurm-apptainer/config.yaml` depending on your HPC and
    workload manager configuration):

    ```bash
    cat << "EOF" > run.sh
    #!/bin/bash
    mkdir -p logs/cluster_log
    snakemake \
        --profile="../profiles/slurm-apptainer" \
        --configfile="config.yaml"
    EOF
    ```

5. Start your workflow run:

    ```bash
    bash run.sh
    ```

6. To find out more information on the expected output files, please visit the
   [output files](./outputs.md) description page.

!!! info "Cluster configuration"

    The [Slurm][slurm] profiles are configured for a cluster that uses the
    quality-of-service (QOS) keyword. If not supported by your Slurm instance,
    you need to remove all the lines with `qos` in
    `profiles/slurm-config.json`. Reach out to your systems administrator for
    further help adjusting the configuration to your specific HPC environment,
    especially when not using the Slurm workload manager.

??? tip "Want to use Conda instead?"

    Change the argument to `--profile` from `local-apptainer` or
    `slurm-apptainer` to `local-conda` or `slurm-conda` instead.

## How to fetch sequencing samples from the Sequence Read Archive?

An independent Snakemake workflow `workflow/rules/sra_download.smk` is included
for the download of sequencing libraries from the [Sequence Read Archive][sra]
and conversion into FASTQ.

The workflow expects the following parameters in the configuration file:

- `samples`: A sample table in TSV format with column *sample*, containing
  *SRR, ERR or DRR identifiers*, as in [this example][sample-table-sra].
- `outdir`: An output directory.
- `samples_out`: The path to the output sample table containing the paths to
  the corresponding outputs FASTQ files.
- `cluster_log_dir`: The directory in which cluster logs are to be stored.

You can use the following call within an activated `zarp` Conda environment to
try an example run of this workflow:

```bash
snakemake \
  --snakefile="workflow/rules/sra_download.smk" \
  --profile="profiles/local-conda" \
  --config \
    samples="tests/input_files/sra_samples.tsv" \
    outdir="results/sra_downloads" \
    samples_out="results/sra_downloads/sra_samples.out.tsv" \
    log_dir="logs" \
    cluster_log_dir="logs/cluster_log"
```

??? question "Want to use Apptainer instead?"

    Change the argument to `--profile` from `local-conda` to `local-apptainer`
    to execute the workflow steps within Apptainer containers.

After successful execution, `results/sra_downloads/sra_samples.out.tsv` should
contain the following:

```tsv
sample  fq1     f2
SRR18552868     results/sra_downloads/compress/SRR18552868/SRR18552868.fastq.gz 
SRR18549672     results/sra_downloads/compress/SRR18549672/SRR18549672_1.fastq.gz       results/sra_downloads/compress/SRR18549672/SRR18549672_2.fastq.gz
ERR2248142      results/sra_downloads/compress/ERR2248142/ERR2248142.fastq.gz 
```

## How to infer sample metadata?

An independent Snakemake workflow `workflow/rules/htsinfer.smk` is available
that populates the sample table required by _ZARP_ with the sample-specific
parameters `seqmode`, `f1_3p`, `f2_3p`, `organism`, `libtype` and `index_size`.
Those parameters are inferred from the provided `fastq.gz` files by
[HTSinfer][htsinfer].

!!! note "Temporary directory"

    The workflow uses the implicit temporary directory from Snakemake, which
    is called with `[resources.tmpdir]`.

The workflow expects the following configuration parameters:

- `samples`: A sample table in TSV format with columns *sample* (containing
  sample identifiers, and *fq1* and *fq2* ()containing the paths to the input
  FASTQ files). See an example [here][sample-table-htsinfer]. If the table
  contains further _ZARP_ compatible columns (see [workflow
  documentation][zarp-workflow-docs], the values specified in these by the user
  are given priority over HTSinfer's results.
- `outdir`: An output directory.
- `samples_out`: The path to the output sample table containing the parameters
  inferred by HTSinfer.
- `records`: The number of sequence records to consider for inference. Set to
  100'000 by default.

You can use the following call within an activated `zarp` Conda environment to
try an example run of this workflow:

```bash
cd tests/test_htsinfer_workflow
snakemake \
  --snakefile="../../workflow/rules/htsinfer.smk" \
  --restart-times=0 \
  --profile="../../profiles/local-conda" \
  --config \
    outdir="results" \
    samples="../input_files/htsinfer_samples.tsv" \
    samples_out="samples_htsinfer.tsv" \
    log_dir="logs" \
    cluster_log_dir="logs/cluster_log" \
  --notemp \
  --keep-incomplete
```

??? question "Want to use Apptainer instead?"

    Change the argument to `--profile` from `local-conda` to `local-apptainer`
    to execute the workflow steps within Apptainer containers.

However, this call will exit with an error, as not all parameters can be
inferred from the example files. The argument `--keep-incomplete` makes sure
the `samples_htsinfer.tsv` file can nevertheless be inspected.

After successful execution - if all parameters could be either inferred or were
specified by the user - `[OUTDIR]/[SAMPLES_OUT]` should contain a populated
table with parameters `seqmode`, `f1_3p`, `f2_3p`, `organism`, `libtype` and
`index_size`.

## How to use Docker?

!!! warning "Experimental feature!"

_ZARP_ is optimised for Linux users as all packages are available via Conda or
Apptainer. Execution on other systems, e.g., Mac OS X, is tricky, even more so
due to the current transition from Intel to ARM processors (M series).
Nevertheless, we built a Docker image that may be used to try to run _ZARP_ in
such environments (tested on Linux and Mac OS).

1. Install Docker following the [official instructions][docker-docs].

2. Pull the Docker image that contains the necessary dependencies:
   ```sh
   docker pull zavolab/zarp:1.0.0-rc.1
   ```

3. Create a directory (e.g. `data/`) and store all the files required for a
   run:
    - The genome sequence FASTA file
    - The annotation GTF file
    - The FASTQ files of your experiments
    - The `rule_config.yaml` for the parameters
    - The `samples.tsv` containing the metadata of your samples
    - The `config.yaml` file with parameters

    Below you can find an example confirguation file. Note how it points to
    files in the `data` directory:

    ```yaml
    ---
    # Required fields
    samples: "data/samples_docker.tsv"
    output_dir: "data/results"
    log_dir: "data/logs"
    cluster_log_dir: "data/logs/cluster"
    kallisto_indexes: "data/results/kallisto_indexes"
    salmon_indexes: "data/results/salmon_indexes"
    star_indexes: "data/results/star_indexes"
    alfa_indexes: "data/results/alfa_indexes"
    # Optional fields
    rule_config: "data/rule_config.yaml"
    report_description: "No description provided by user"
    report_logo: "../../images/logo.128px.png"
    report_url: "https://zavolan.biozentrum.unibas.ch/"
    author_name: "NA"
    author_email: "NA"
    ...
    ```

4. Execute _ZARP_:
    ```sh
    docker run \
        --platform linux/x86_64 \
        --mount type=bind,source=$PWD/data,target=/data \
        zavolab/zarp:1.0.0-rc.1 \
        snakemake -p \
        --snakefile /workflow/Snakefile \
        --configfile data/config.yaml \
        --cores 4 \
        --use-conda \
        --verbose
    ```

    The command runs the Docker container `zavolab/zarp:1.0.0-rc.1` that we
    have pulled. It executes it, just as it would be done on a Linux platform
    (`--platform linux/x86_64`). We use the `--mount` option to bind the local
    `data/` directory that contains the input files to the `/data/` directory
    inside the container. The workflow definition file is stored at
    `/workflow/Snakefile` inside the container.

    Once _ZARP_ has completed, the results will be available in the
    `data/results/` directory.
