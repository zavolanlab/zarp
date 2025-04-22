# Parameterization

## Adjusting tool parameters

!!! warning "Experimental feature!"

To increase usability and ensure _ZARP_ produces reasonable results for the
vast majority of RNA-Seq samples that it was built, we have consciously
limited the set of tool configuration parameters that the workflow exposes
through Snakemake's `config.yaml`.

However, recognizing that power users would surely like to tweak tool
parameters (e.g., to allow _ZARP_ to be run against non-standard RNA-Seq
protocols), we provide an additional, custom config file `rule_config.yaml`,
which enables the user to modify (almost) any tool parameter. The values
specified in this configuration override the default ones, unless they are
deemed essential to the correct "wiring" (in which case they are considered
"immutable" and will not be changed).

Modifying `rule_config.yaml` thus allows users to fundamentally alter _ZARP_'s
behavior, while keeping its general wiring, all _without_ having to make
changes in the workflow definition itself.

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

You can find the path to the `rule_config.yaml` file as a parameter in the
standard `config.yaml` file.

```bash
rule_config: "../input_files/rule_config.yaml"
```

## Upgrading tool versions

As described elsewhere, _ZARP_'s rules can be executed either with
[Conda][conda] or with [Apptainer][apptainer]. The majority of the tools that
we use are hosted by the [Bioconda][bioconda] and
[BioContainers][biocontainers] registries, for Conda environments and container
images, respectively.

Incase you want to upgrade one of the tools to a later version, all you need
to do is to update the conda environment "recipe" file in `workflow/envs/` and
then modify the corresponding container directive in the appropriate Snakemake
workflow definition files (either in `workflow/Snakefile` or in one of the
imported "subworkflows" in `workflow/rules`).

??? warning "We do not recommend downgrading tool versions!"

    For security and compatibility reasons, we strongly recommend only to
    _upgrade_ tool versions!

For example, let's say you want to update `cutadapt` to the latest version:

1. Find the version you want to use by searching for the correspondign
   tool/package name on the [Anaconda website][anaconda]. For `cutadapt`,
   the page for the specific Bioconda package is available
   [here][bioconda-cutadapt]. At the time of writing, the latest version
   is `4.9`.

    <div align="center">
        <img width="80%" src=../images/bioconda_cutadapt.png>
    </div>

2. Find the appropriate Conda enviroment recipe file in `workflow/envs/`. In
   our example, it is `workflow/envs/cutadapt`, and it contains the following:

    ```yaml
    ---
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - cutadapt=4.6
    ...
    ```

    Simply replace `4.6` with the `4.9` you identified in step above.

3. Now find the corresponding BioContainers container image from the [quay.io
   website][quay]. Select the tags and note down one of the available versions
   corresponding to `cutadapt 4.9` (note that they are available for multiple
   Python versions; generally pick the most recent one you are comfortable
   with). Copy the _entire_ tag name, not just the `4.9` part!

    <div align="center">
        <img width="80%" src=../images/biocontainers_cutadapt.png>
    </div>

4. Finally, determine all the places where `cutadapt` is used in the workflow,
   by inspecting `workflow/Snakefile` and all of the `.smk` files in
   `workflow/rules/`. In the case of `cutadapt`, multiple rules use the tool,
   each with the following `container` and `conda` directives.

    ```
    container:
        "docker://quay.io/biocontainers/cutadapt:4.6--py310h4b81fae_1"
    conda:
        os.path.join(workflow.basedir, "envs", "cutadapt.yaml")
    ```

    For each of these, replace the final part of the `container` directive
    (here: `4.6--py310h4b81fae_1`) with the tag name you copied in the step
    above.

That's it.

## Adding new tools

You may have found or developed a package that you would like to include in
_ZARP_, so that it is always executed whenever you start a run. In this case,
we highly recommend packaging the tool properly and making it available on
Bioconda (see [here][bioconda-contributing] for instructions). The advantage
of this approach is that it allows you to share your tool with the broader
scientific community, ensuring that it is easily accessible and installable by
others. Once your package is available on Bioconda, BioContainers will
automatically build a Docker image for your package (available via
[quay.io][quay]) - which you Apptainer will be able to use. With your Bioconda
package and container image available, you can easily add additional rules to
the Snakemake definition files, inside your own copy of the _ZARP_ repository.

Or perhaps you are convinced that every _ZARP_ user should always run the tool
or tools you have added? In that case, create a pull request against the
upstream/original _ZARP_ repository. We will evaluate your request and - if we
like it and tests pass - merge it.

Don't know how to do that? Just write us a brief [email][contact]! :relaxed:

## Managing tool resources

_ZARP_ has been tested with many samples, and we set default parameters to
ensure optimal performance of the tools across a broad range of samples (e.g.,
a tool can run on multiple threads). With respect to setting the maximum
memory usage per rule, Snakemake provides a variable `mem_mb` in the directive
`resources`. But given that the memory usage is strongly sample-dependent
(sample size, souce organism), in _ZARP_ we use _dynamic memory allocation_.

We estimate the required memory based on the size of the input file using
scaling factors that we have obtained empirically from the analysis of many
samples. Should that initial estimate lead to a rule failing, it is re-run
with an increased memory allocation. This process is repeated up to three
times.

This is what the _dynamic memory allocation_ looks like in the workflow
definition file:

```
resources:
    mem_mb=lambda wildcards, attempt: 4096 * attempt,
```

In some cases the workflow might still fail. This means that you would need to
manually alter the `mem_mb` used. The best solution for that is to increase the
orginal memory used. In the above example, increase it from `4096` MB to
a higher value, based on your expectations.

??? note "Globally capping resource consumption"

    Note that Snakemake provides configuration parameters to globally limit
    resource usage (both for memory and the number of CPUs/threads). If set,
    these will take precedence over the individual rule settings. This is
    useful if you are working on a laptop, or another machine with very
    limited resources. However, if your resources are lower than what the
    most resource-hungry tool integrated in _ZARP_ requires as a minimum for
    a given run, you will not be able to complete that run.

When you submit a workflow to an HPC cluster, additional parameters can be
configured through [Snakemake profiles][snakemake-profiles]. In the
`profiles/` directory, yo can a few basic options that have proven useful for
us. For example, the `slurm-conda` profile is available, and, as the name
suggests, it submits jobs to a [Slurm][slurm]-managed HPC cluster while using
[Conda][conda] to manage the dependencies for each workflow rule.

Using profiles, you can also set default resources that are applied to all
jobs, unless explicitly specified. For example, if you want to set the default
memory used for all rules, add or modify the following line in the profile
configuration:

```
default-resources: mem_mb=1024
```

??? question "Where to configure runtime limits?"

    Some parameters, including runtime limits, can only be set in your
    workload manager configuration. In the case of Slurm, this is in
    `slurm-config.json`.
