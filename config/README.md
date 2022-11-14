# Running the workflow on your own samples

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

    When running the pipeline with *conda* you should use `local-conda` and
    `slurm-conda` profiles instead.

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
