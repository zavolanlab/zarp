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
    rm -rf snakemake_report.html
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
    --snakefile="../../workflow/Snakefile" \
    --profile="../../profiles/local-singularity" \
    --configfile="../input_files/config.yaml"

# Create a Snakemake report after the workflow execution
snakemake \
    --snakefile="../../workflow/Snakefile" \
    --configfile="../input_files/config.yaml" \
    --report="snakemake_report.html"

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output_temp_flag.md5"

# Checksum file generated with
# find results/ \
#    -type f \
#    -name \*\.gz \
#    -exec gunzip '{}' \;
# find results/ \
#    -type f \
#    -name \*\.zip \
#    -exec sh -c 'unzip -o {} -d $(dirname {})' \;
# md5sum $(cat expected_output.files) > expected_output_temp_flag.md5

# Check whether STAR produces expected alignments
# STAR alignments need to be fully within ground truth alignments for tests to pass; not checking 
# vice versa because processing might cut off parts of reads (if testing STAR directly, add '-f 1' 
# as additional option)
echo "Verifying STAR output"
result=$(bedtools intersect -F 1 -v -bed \
    -a ../input_files/synthetic.mate_1.bed \
    -b results/samples/synthetic_10_reads_mate_1_synthetic_10_reads_mate_1/map_genome/synthetic_10_reads_mate_1_synthetic_10_reads_mate_1.se.Aligned.sortedByCoord.out.bam \
    | wc -l)
if [ $result != "0" ]; then
    echo "Alignments for mate 1 reads are not consistent with ground truth"
    exit 1
fi
result=$(bedtools intersect -F 1 -v -bed \
    -a <(cat ../input_files/synthetic.mate_1.bed ../input_files/synthetic.mate_2.bed) \
    -b results/samples/synthetic_10_reads_paired_synthetic_10_reads_paired/map_genome/synthetic_10_reads_paired_synthetic_10_reads_paired.pe.Aligned.sortedByCoord.out.bam \
    | wc -l)
if [ $result != "0" ]; then
    echo "Alignments for mate 1 reads are not consistent with ground truth"
    exit 1
fi

# Check whether Salmon assigns reads to expected genes
echo "Verifying Salmon output"
diff \
    <(cat results/samples/synthetic_10_reads_mate_1_synthetic_10_reads_mate_1/synthetic_10_reads_mate_1_synthetic_10_reads_mate_1.salmon.se/quant.genes.sf | cut -f1,5 | tail -n +2 | sort -k1,1) \
    <(cat ../input_files/synthetic.mate_1.bed | cut -f7 | sort | uniq -c | sort -k2nr | awk '{printf($2"\t"$1"\n")}')
diff \
    <(cat results/samples/synthetic_10_reads_paired_synthetic_10_reads_paired/synthetic_10_reads_paired_synthetic_10_reads_paired.salmon.pe/quant.genes.sf | cut -f1,5 | tail -n +2 | sort -k1,1) \
    <(cat ../input_files/synthetic.mate_1.bed | cut -f7 | sort | uniq -c | sort -k2nr | awk '{printf($2"\t"$1"\n")}')

