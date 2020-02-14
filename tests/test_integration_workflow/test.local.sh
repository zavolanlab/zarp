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
mkdir -p logs/local_log

# Run tests
snakemake \
    --snakefile="../../snakemake/Snakefile" \
    --configfile="config.yaml" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args="--bind ${PWD}"

# Check md5 sum of some output files
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

# Check whether STAR produces expected alignments
# STAR alignments need to be fully within ground truth alignments for tests to pass; not checking 
# vice versa because processing might cut off parts of reads (if testing STAR directly, add '-f 1' 
# as additional option)
echo "Verifying STAR output"
result=$(bedtools intersect -F 1 -v -bed \
    -a input_files/synthetic.mate_1.bed \
    -b results/single_end/synthetic_10_reads_mate_1/map_genome/synthetic_10_reads_mate_1_Aligned.sortedByCoord.out.bam \
    | wc -l)
if [ $result != "0" ]; then
    echo "Alignments for mate 1 reads are not consistent with ground truth"
    exit 1
fi
result=$(bedtools intersect -F 1 -v -bed \
    -a input_files/synthetic.mate_2.bed \
    -b results/single_end/synthetic_10_reads_mate_2/map_genome/synthetic_10_reads_mate_2_Aligned.sortedByCoord.out.bam \
    | wc -l)
if [ $result != "0" ]; then
    echo "Alignments for mate 1 reads are not consistent with ground truth"
    exit 1
fi
result=$(bedtools intersect -F 1 -v -bed \
    -a <(cat input_files/synthetic.mate_1.bed input_files/synthetic.mate_2.bed) \
    -b results/paired_end/synthetic_10_reads_paired/map_genome/synthetic_10_reads_paired_Aligned.sortedByCoord.out.bam \
    | wc -l)
if [ $result != "0" ]; then
    echo "Alignments for mate 1 reads are not consistent with ground truth"
    exit 1
fi

# Check whether Salmon assigns reads to expected genes
echo "Verifying Salmon output"
diff \
    <(cat results/single_end/synthetic_10_reads_mate_1/salmon_quant/quant.genes.sf | cut -f1,5 | tail -n +2 | sort -k1,1) \
    <(cat input_files/synthetic.mate_1.bed | cut -f7 | sort | uniq -c | sort -k2nr | awk '{printf($2"\t"$1"\n")}')
diff \
    <(cat results/single_end/synthetic_10_reads_mate_2/salmon_quant/quant.genes.sf | cut -f1,5 | tail -n +2 | sort -k1,1) \
    <(cat input_files/synthetic.mate_2.bed | cut -f7 | sort | uniq -c | sort -k2nr | awk '{printf($2"\t"$1"\n")}')
diff \
    <(cat results/paired_end/synthetic_10_reads_paired/salmon_quant/quant.genes.sf | cut -f1,5 | tail -n +2 | sort -k1,1) \
    <(cat input_files/synthetic.mate_1.bed | cut -f7 | sort | uniq -c | sort -k2nr | awk '{printf($2"\t"$1"\n")}')

