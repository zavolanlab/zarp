# Parameterization

## Parameter adjustment for rules

ZARP runs on default parameters that were chosen based on the majority of samples being analyzed. To simplify adding non-default parameters of the tools, we created a config file, which enables the user to use any of the extra parameters offered by any of the tools. These options override the default ones, unless they are essential to the correct use of the rule in which case they are considered "immutable" and will not be changed. This extra `rule_config.yaml` file enables the user to easily manipulate the parameters, without having to deal with making changes in the workflow.

An example of a `rule_config.yaml` is shown below:

```bash
remove_adapters_cutadapt:
    # Search for all the given adapter sequences repeatedly, either until no
    # adapter match was found or until n rounds have been performed (default 1,
    # ZARP recommends 2)
    -n: '2'
    # Discard processed reads that are shorter than m; note that cutadapt uses
    # a default value of m=0, causing reads without any nucleotides remaining
    # after processing to be retained; as "empty reads" will cause errors in
    # downstream applications in ZARP, we have changed the default to m=1,
    # meaning that only read fragments of at least 1 nt will be retained after
    # processing. The default will be overridden by the value specified here,
    # but for the reason stated above, we strongly recommend NOT to set m=0;
    # cf. https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
    -m: '10'
```

You can find the path to the `rule_config.yaml` file as a parameter in the standard `config.yaml` file.

```bash
rule_config: "../input_files/rule_config.yaml"
```

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
    
2. Find the version you want to use by searching the package in the anaconda website. The specific package is available [here](https://anaconda.org/bioconda/cutadapt). The latest version at the moment is 4.9.

    <div align="center">
        <img width="80%" src=../images/bioconda_cutadapt.png>
    </div>

3. Find the corresponding version from biocontainers. You can do that by searching something like "biocontainers cutadapt" in the [quay.io](https://quay.io/) website. Select the tags and use one of the available versions.

    <div align="center">
        <img width="80%" src=../images/biocontainers_cutadapt.png>
    </div>

4. You can replace the dependencies with the new versions.

## Expand zarp with a new package

In other cases, you might have developed a package and want to expand ZARP with it. In that case, we highly recommend packaging the tool properly and making it available on Bioconda. You can find the guidelines [here](https://bioconda.github.io/contributor/index.html). This is how we also publish custom tools developed in the lab. The advantage of this approach is that it allows you to share your tool with the broader scientific community, ensuring that it is easily accessible and installable by others. Once your package is available on Bioconda, you will also get a Docker container that contains your package. This is automatically built by the Biocontainers team and will become available on [quay.io](https://quay.io/) as shown in the previous section.

## Update the tool resources

ZARP has been tested with many samples and we provided default parameters to allow optimal performance of the tools (e.g., a tool can run on multiple threads). Regarding the required maximum memory usage required, there is a field in the snakemake rule called resources where the maximum memory per rule can be specified using the variable mem_mb. In the case of ZARP we go one step further and use dynamic resources, which means we estimate the required memory based on the size of the input file using scaling factors that we have obtained from analysis of multiple samples. An extra dynamic modification is that if the rule fails, we allow three reruns of the rule during which we multiply the provided memory by the attempt. That means that if the rule fails, in the rerun the memory will be doubled or tripled. This is done as following:

```
resources:
    mem_mb=lambda wildcards, attempt: 4096 * attempt,
```

In some cases the pipeline might still fail. This means that you would need to alter the mem_mb used. The best solution for that is to increase the orginal memory used. In the above example increase it from `4096` MB to something higher. If the user limits the available overall memory when executing snakemake through the resources option, this overrides the per rule specifications. All of this depends on the user system having the required resources. 

Similarly in each rule the number of cores is specified via the `threads` parameter, however there is a global snakemake specification of cores, which overrides the per rule specifications and might be limiting the efficiency of the workflow.

When you submit a workflow to a HPC cluster additional parameters can be customized through `profiles`. By default under the `profiles` directory you can find different options. For example the `slurm-conda` profile is available and as the name suggests submits jobs to a slurm cluster and uses conda for the package dependencies. If you want to increase the default memory used for all the jobs, this can easily happen by altering the following line:
```
default-resources: mem_mb=1024
```
Other parameters (e.g., time) can be customized in the slurm `slurm-config.json` file.