# set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
--cores 4 \
-p \
--rerun-incomplete \
--use-singularity \
--singularity-args "--bind ${PWD}/../tests/input_files"
