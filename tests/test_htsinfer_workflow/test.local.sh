#!/bin/bash

# This script is currently exiting with non-zero status.
# This is expected behaviour though, as several parameters can't be inferred from the test files.
# Snakemake called with --keep-incomplete in order to keep the created samples table for inspection.

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf .cache/
    rm -rf .config/
    rm -rf .fontconfig/
    rm -rf .java/
    rm -rf .snakemake/
    rm -rf logs/
    rm -rf results/
    rm -rf Log.out
    cd $user_dir
    echo "Exit status: $rc"
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

# Run tests
snakemake \
    --snakefile="../../workflow/rules/htsinfer.smk" \
    --restart-times=0 \
    --profile="../../profiles/local-singularity" \
    --config outdir="results" \
             samples="../input_files/htsinfer_samples.tsv" \
             samples_out="samples_htsinfer.tsv" \
    --notemp \
    --keep-incomplete

