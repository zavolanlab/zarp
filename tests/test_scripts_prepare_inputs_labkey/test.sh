#!/bin/bash

# Scripts requires environment variables 'LABKEY_HOST', 'LABKEY_USER' and
# 'LABKEY_PASS' to be set with the appropriate values

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf ${HOME}/.netrc
    rm -rf .snakemake/
    rm -rf config.yaml
    rm -rf samples.tsv.labkey
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
python "../../scripts/prepare_inputs.py" \
    --labkey-domain="${LABKEY_HOST}" \
    --labkey-path="/Zavolan Group/TEST_LABKEY" \
    --input-to-output-mapping="../../scripts/prepare_inputs.dict.tsv" \
    --resources-dir="../input_files" \
    --output-table="samples.tsv" \
    --config-file="config.yaml" \
    --multimappers='10' \
    --logo="../../images/logo.128px.png" \
    --debug \
    "RNA_Seq_data_template_raw"

# Check if dry run completes
snakemake \
    --snakefile="../../Snakefile" \
    --configfile="config.yaml" \
    --dryrun \
    --verbose

#md5sum --check "expected_output.md5"
# MD5 sums obtained with command:
# md5sum config.yaml samples.tsv > expected_output.md5
md5sum config.yaml samples.tsv
