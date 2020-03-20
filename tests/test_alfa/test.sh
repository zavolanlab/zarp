#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf .java/
    rm -rf .snakemake/
    rm -rf logs/
    rm -rfv results/alfa_indexes/
    rm -rfv results/*/*/ALFA/
    rm -rfv results/*/*/STAR_coverage/
    rm -rfv results/ALFA/
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
    --snakefile="../../Snakefile" \
    --configfile="../input_files/config_alfa.yaml" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../input_files,${PWD}/../../images" \
    --verbose \
    results/ALFA/ALFA_plots.concat.png

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"
