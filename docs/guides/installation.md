# Installation

On this page, you will find out how to install _ZARP_ on your system.

## Requirements

Installation requires the following:

- Linux (tested with Ubuntu 20.04; macOS has not been tested yet)
- [Conda][conda] (tested with `Conda 22.11.1`)
- [Mamba][mamba] (tested with `Mamba 1.3.0`)
- [Singularity][singularity] (tested with `Singularity 3.8.6`; not required
  if you have root permissions on the machine you would like to install _ZARP_
  on; in that case, see [below](#2-set-up-conda-environment))

> Other versions, especially older ones, are not guaranteed to work.

## Installation steps

### 1. Clone ZARP

Clone the [ZARP workflow repository][zarp] with:

```sh
git clone https://github.com/zavolanlab/zarp.git
```

### 2. Set up Conda environment

To check if you already have conda installed in your system type:

```sh
conda --version
```

Example:
```sh
conda --version
conda 22.11.1
```
If it is not installed, you will see a <code style="color: red;">command not found error.</code>

#### Conda installation

There are different ways to install Conda. We recommend to use [Miniconda][miniconda] to install it. Please refer to the official documentation for the latest installation. Alternatively you can follow these steps:

**1. Download the Miniconda installer:**

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh
```
**2. Run the installer:**
```sh
bash Miniconda3-4.7.12-Linux-x86_64.sh
```
**3. Follow the prompts** to complete the installation

**4. Initialize Conda:**
```sh
source ~/.bashrc
```
After installation, you can verify it again using ```conda --version.```

**5. Update Conda to a specific version (e.g., 22.11.1):**
```sh
conda install conda=22.11.1
```
>This update includes a step to install a specific version of Conda, ensuring that users have a version tested to be compatible with ZARP.

#### Conda installation if you already have Conda and do NOT want to change its version

If you already have a specific conda version on your system which is not compatible with ZARP and do not want to change it, no worries. You can have more than two conda versions:

**1. Download the installer for the second version of Conda:**

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.7.12-Linux-x86_64.sh
```
**2.Install it in a different directory:**
```sh
bash Miniconda3-py39_4.7.12-Linux-x86_64.sh -b -p $HOME/miniconda3_py39
```
**3.Source the appropriate conda.sh file to switch between versions:**

To use the first version:
```sh
source $HOME/miniconda3/etc/profile.d/conda.sh
```
To use the second version:
```sh
source $HOME/miniconda3_py39/etc/profile.d/conda.sh
```

**4. Update Conda to a specific version (e.g., 22.11.1):**
```sh
conda install conda=22.11.1
```

### 3. Set up mamba
Given that Miniconda has been installed and is available in the current shell the first dependency for ZARP is the [Mamba][mamba] package manager, which needs to be installed in the base conda environment with:
```sh
conda install mamba=1.3.0 -n base -c conda-forge
```

### 4. Set up your ZARP environment

For improved reproducibility and reusability of the workflow, each individual step of the workflow runs either in its own Singularity container or in its own Conda virtual environemnt. As a consequence, running this workflow has very few individual dependencies. The container execution requires Singularity to be installed on the system where the workflow is executed. As the functional installation of Singularity requires root privileges, and Conda currently only provides Singularity for Linux architectures, the installation instructions are slightly different depending on your system/setup:

**For most users**

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

**As root user on Linux**

If you have a Linux machine, as well as root privileges, (e.g., if you plan to
run the workflow on your own computer), you can execute the following command
to include Singularity in the Conda environment:

```bash
mamba env update -f install/environment.root.yml
```

### 5. Activate ZARP environment

Activate the Conda environment with:

```bash
conda activate zarp
```

### 6. Optional installation steps

#### Install test dependencies

Most tests have additional dependencies. If you are planning to run tests, you
will need to install these by executing the following command _in your active
Conda environment_:

```bash
mamba env update -f install/environment.dev.yml
```

#### Run installation tests

We have prepared several tests to check the integrity of the workflow and its
components. These can be found in subdirectories of the `tests/` directory. 
The most critical of these tests enable you to execute the entire workflow on a 
set of small example input files. Note that for this and other tests to complete
successfully, [additional dependencies](#installing-non-essential-dependencies) 
need to be installed. 
Execute one of the following commands to run the test workflow 
on your local machine:


##### Test workflow on local machine with **Singularity**

```bash
bash tests/test_integration_workflow/test.local.sh
```

##### Test workflow on local machine with **Conda**

```bash
bash tests/test_integration_workflow_with_conda/test.local.sh
```
Execute one of the following commands to run the test workflow 
on a [Slurm][slurm]-managed high-performance computing (HPC) cluster:

##### Test workflow with **Singularity**

```bash
bash tests/test_integration_workflow/test.slurm.sh
```

##### Test workflow with **Conda**

```bash
bash tests/test_integration_workflow_with_conda/test.slurm.sh
```

> **NOTE:** Depending on the configuration of your Slurm installation you may
> need to adapt file `slurm-config.json` (located directly under `profiles`
> directory) and the arguments to options `--cores` and `--jobs`
> in the file `config.yaml` of a respective profile.
> Consult the manual of your workload manager as well as the section of the
> Snakemake manual dealing with [profiles].

