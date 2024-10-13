# Contribute

Open source contributors are always welcome, for [_ZARP_][zarp], [_ZARP-cli_][zarp-cli] or any other of the [Zavolab projects][zavolab-gh]. Simply reach out by [email][contact] to schedule an onboarding call.

## Update the version of a tool

Each step (rule) of ZARP can be executed either with conda, or with apptainer (singularity). The majority of the tools that we use are hosted by the [Bioconda](https://bioconda.github.io/) and [BioContainers](https://biocontainers.pro/registry) registries. So in case you want to update one of the tools into a later version you need to update both the conda yaml file of the rule and the corresponding container. For example let's say you want to update cutadapt to the latest version.

1. Determine the place where `cutadapt` is used. Multiple rules use the following dependencies:
    ```
    container:
        "docker://quay.io/biocontainers/cutadapt:4.6--py310h4b81fae_1"
    conda:
        os.path.join(workflow.basedir, "envs", "cutadapt.yaml")
    ```

    The `cutadapt.yaml` looks like the following:
    ```yaml
    ---
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - cutadapt=4.6
    ...
    ```
    
2. Find the version you want to use by searching the package in the anaconda website. The specific package is available [here](https://anaconda.org/bioconda/cutadapt).

<div align="center">
    <img width="80%" src=../images/bioconda_cutadapt.png>
</div>

The latest version at the moment is 4.9.

3. Find the corresponding version from biocontainers. You can do that by searching something like "biocontainers cutadapt" in the [quay.io](https://quay.io/) website. Select the tags and use one of the available versions.

<div align="center">
    <img width="80%" src=../images/biocontainers_cutadapt.png>
</div>

4. You can replace the dependencies with the new versions.