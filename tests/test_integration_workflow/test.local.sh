#!/bin/bash

# Tear down test environment
trap 'rm -rf logs/ results/ .snakemake/ && cd $user_dir' EXIT  # quotes command is exected after script exits, regardless of exit status

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir
mkdir -p logs/local_log

# Run tests
snakemake \
    --snakefile="../../snakemake/Snakefile" \
    --configfile="config.yaml" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args "--bind ${PWD}"
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
md5sum --check "expected_output.md5"

# Checksum file generated with
# find results/ \
#     -type f \
#     -name \*\.gz \
#     -exec gunzip '{}' \;
# find results/ \
#     -type f \
#     -regextype posix-egrep \
#     -regex ".*\.(fastq|html)$" \
#     -exec md5sum '{}' \; \
#     > expected_output.md5

