#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf .cache/
    rm -rf .config/
    rm -rf .fontconfig/
    rm -rf .java/
    rm -rf .snakemake/
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

mkdir -p data
cp ../../tests/input_files/homo_sapiens/genome.fa data/genome.fa
cp ../../tests/input_files/homo_sapiens/annotation.gtf data/annotation.gtf
cp -r ../../tests/input_files/project1 data/project1
cp -r ../../tests/input_files/project2 data/project2
cp -r ../../tests/input_files/config_docker.yaml data/config_docker.yaml
cp ../../tests/input_files/rule_config.yaml data/rule_config.yaml
cp ../../tests/input_files/samples_docker.tsv data/samples_docker.tsv

# Pull the zarp container
docker pull zavolab/zarp:1.0.0-rc.1

# Run tests with Docker
docker run \
    --platform linux/x86_64 \
    --mount type=bind,source=$script_dir/data,target=/data \
    zavolab/zarp:1.0.0-rc.1 \
    snakemake \
    -p \
    --snakefile /workflow/Snakefile \
    --configfile data/config_docker.yaml \
    --cores 4 --use-conda --verbose
