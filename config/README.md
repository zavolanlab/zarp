# Dependencies installation

Running this workflow has very few individual dependencies. 

```bash
mamba env create -f install/environment.yml
```

# Activate environment

Activate the Conda environment with:

```bash
conda activate zarp
```

# Run the workflow on your own samples

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

    - [samples.tsv](../tests/input_files/samples.tsv)
    - [config.yaml](../tests/input_files/config.yaml)


# More execution options
For more execution options, like sample fetching from SRA, sample features inferense, 
, cluster execution and many more, visit ZARP [../README.md] and [../pipeline-documentation]

