#!/bin/bash

# Tear down test environment
trap 'rm -rf logs/ results/ .snakemake/ && cd $user_dir' EXIT  # quoted command is exected after script exits, regardless of exit status

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir
mkdir -p logs/cluster_log
mkdir -p logs/local_log

# Run tests
snakemake \
    --snakefile="../../snakemake/Snakefile" \
    --configfile="config.yaml" \
    --cluster-config="cluster.json" \
    --cluster="sbatch --cpus-per-task={cluster.threads} --mem={cluster.mem} --qos={cluster.queue} --time={cluster.time} --job-name={cluster.name} -o {cluster.out} -p scicore" \
    --cores=256 \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args="--bind ${PWD}"
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"

# Checksum file generated with
# find results/ \
#     -type f \
#     -name \*\.gz \
#     -exec gunzip '{}' \;
# find results/ \
#     -type f \
#     -name \*\.zip \
#     -exec sh -c 'unzip -o {} -d $(dirname {})' \;
# md5sum $(cat expected_output.files) > expected_output.md5

