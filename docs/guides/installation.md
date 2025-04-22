# Installation

On this page, you will find out how to install _ZARP_ on your system.

## Requirements

For improved reproducibility and reusability, each individual workflow step
runs either inside a dedicated [Conda][conda] environment or an
([Apptainer][apptainer]) container. As a consequence, running _ZARP_ only has
very few dependencies, that need to be available on your system:

- Linux (tested with `Ubuntu 24.04`)
- [Conda][conda] (tested with `Conda 24.11.3`)
- **Optional:** [Apptainer][apptainer] (tested with `Apptainer 1.3.6`)

A few additional dependencies are installed via Conda as described further
below.

!!! warning "Other versions, especially older ones, are not guaranteed to work!"

??? question "Don't have Linux?"

    Please see the "How to use Docker?" instructions in the [usage
    section](./usage.md).

??? question "How do I install Apptainer?"

    Please follow the [official documentation][apptainer-docs] to install
    Apptainer (formerly Singularity) globally and configure its permissions.

## Installation steps

### 1. Clone _ZARP_

Clone the [_ZARP_ workflow repository][zarp] with:

```sh
git clone https://github.com/zavolanlab/zarp.git
```

### 2. Ensure Conda is available

To check if you already have Conda installed on your system, type:

```sh
conda --version
```

If Conda is available, you should see output similar to this:

```console
conda 24.11.3
```
If it is not installed, you will instead see
<code style="color: red;">command not found: conda</code>.

??? question "What's the best way to install Conda?"

    Conda can be installed in multiple ways. We strongly recommend using the
    [Miniforge][miniforge] distribution, as it is built around the
    community-supplied `conda-forge` channel that is heavily used by _ZARP_.
    It also reduces the risk of accidentally violating Conda's licensing
    restrictions by using the `defaults` channel.

    Please refer to the [official documentation][miniforge] for up-to-date
    installation instructions.

??? tip "My Conda version is not compatible"

    After completing Conda setup, you can install a specific Conda version with
    the following command:

    ```sh
    conda install conda=24.11.3
    ```

??? tip "I do not want to change my Conda version"

    If you already have a specific Conda version on your system that is not
    compatible with _ZARP_, and you do not want to change it, no worries. Just
    indicate a different directory when using the interactive Miniforge
    installer. Then source the appropriate `conda.sh` file to switch between
    versions, e.g.:

    ```sh
    source $HOME/miniconda3/etc/profile.d/conda.sh  # OR
    source $HOME/miniconda3_alt/etc/profile.d/conda.sh
    ```

### 3. Set up your _ZARP_ environment

To install the remaining _ZARP_ dependencies, run:

```bash
conda env create -f install/environment.yml
```

??? tip "Installing development dependencies"

    If you would like to [run tests](#running-installation-tests) or if you
    are planning to contribute to _ZARP_'s development, run the following
    command instead to set up your _ZARP_ Conda environment.

    ```bash
    conda env create -f install/environment.dev.yml
    ```

    This will ensure that all development dependencies are installed as well.

??? tip "You want to run _ZARP_ on an HPC?"

    When running _ZARP_ on a High-Performance Computing cluster, you will need
    to make sure that compatible versions of at least one of Conda (when using
    Snakemake's `--use-conda` option) and Apptainer (when using the
    `--use-apptainer` option) are installed and properly configured on each
    machine of the cluster, including the head node. Reach out to your systems
    administrator to set this up for you.

### 4. Activate the _ZARP_ environment

Activate the Conda environment with:

```bash
conda activate zarp
```

## Running installation tests

We have prepared several tests to check the integrity of the workflow and its
components. These can be found in subdirectories of the `tests/` directory.
The most critical of these tests enable you to execute the entire workflow on a
set of small example input files. Note that for this and other tests to
complete successfully, additional [development
dependencies](#3-set-up-your-zarp-environment) need to be installed.

Execute the commands below to run the test workflow on your local machine or
on your [Slurm][slurm]-managed HPC cluster.

??? tip "Failing tests do not necessarily indicate a problem!"

    Our tests were developed to guard against code regression over time or as
    a result of proposed changes and are therefore very rigorous. In
    particular, even the minutest changes in outputs, even individual pixels
    or metadata in the produced output images, will cause a test suite to
    fail. However, we cannot rule out such minor changes across systems, and
    therefore, it is quite possible that your test runs may fail, even though
    the workflow is properly installed and functional. Check the logs to see
    whether Snakemake completed successfully. If so, it is very likely that
    everything is fine.

### Running tests on your local machine

Use the following command to run each step inside a dedicated **Conda
environment**:

```bash
bash tests/test_integration_workflow_with_conda/test.local.sh
```

Instead, use the following command to run each step inside an **Apptainer
container**:

```bash
bash tests/test_integration_workflow_with_conda/test.slurm.sh
```

### Running tests on your Slurm cluster

Use the following command to run each step inside a dedicated **Conda
environment**:

```bash
bash tests/test_integration_workflow/test.local.sh
```

Instead, use the following command to run each step inside an **Apptainer
container**:

```bash
bash tests/test_integration_workflow/test.slurm.sh
```

??? tip "The Slurm tests are failing for me!"

    Depending on the configuration of your Slurm installation you may need to
    adapt file `profiles/slurm-config.json` and the arguments to options
    `--cores` and `--jobs` in the file `config.yaml` of a respective profile.
    Consult the manual of your workload manager as well as the section of the
    Snakemake manual dealing with [profiles][snakemake-profiles].
