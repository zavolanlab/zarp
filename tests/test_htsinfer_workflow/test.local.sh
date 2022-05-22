#!/bin/bash

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
    --profile="../../profiles/local-singularity-htsinfer" \
    --config outdir="results" samples="../input_files/htsinfer_samples.tsv" samples_out="samples_htsinfer.tsv" \
    --notemp

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"

