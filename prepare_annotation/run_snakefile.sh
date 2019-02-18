# set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
--cluster-config cluster.json \
--cluster "sbatch --cpus-per-task={cluster.threads} --mem={cluster.mem} --qos={cluster.queue} --time={cluster.time} --job-name={cluster.name} -o {cluster.out} -p scicore" \
--cores 256 \
-p \
--rerun-incomplete \
--use-singularity \
--singularity-args "--bind ${PWD}"
