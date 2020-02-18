#!/bin/bash

# Scripts requires environment variables 'LABKEY_HOST', 'LABKEY_USER' and
# 'LABKEY_PASS' to be set with the appropriate values

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf ${HOME}/.netrc
    rm -rf .snakemake/
    rm -rf config.yaml
    rm -rf input_table.tsv
    rm -rf samples.tsv
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

cat << EOF | ( umask 0377; cat >> ${HOME}/.netrc; )
machine ${LABKEY_HOST}
login ${LABKEY_USER}
password ${LABKEY_PASS}
EOF

# Run tests
python "../../scripts/labkey_to_snakemake.py" \
    --input-dict="../../scripts/labkey_to_snakemake.dict.tsv" \
    --config-file="config.yaml" \
    --samples-table="samples.tsv" \
    --multimappers='10' \
    --remote \
    --project-name "TEST_LABKEY" \
    --table-name "RNA_Seq_data_template" \
    "../input_files"

# Check if dry run completes
snakemake \
    --snakefile="../../Snakefile" \
    --configfile="config.yaml" \
    --dryrun \
    --verbose

md5sum --check "expected_output.md5"
# MD5 sums obtained with command:
# md5sum config.yaml samples.tsv > expected_output.md5
