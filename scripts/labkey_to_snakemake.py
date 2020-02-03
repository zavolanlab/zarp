## -----------------------------------------------------------------------------
# Author : Katsantoni Maria, Christina Herrmann
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This script is part of the Zavolan lab quantification pipeline, which is used
# for analysing RNA-seq data. The table is provided by labkey and is a csv file.
# If the user provides their own table the table should contain the following 
# columns:
# -----------------------------------------------------------------------------

import sys
import gzip
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO
from csv import writer
from pathlib import Path

# ----------------------------------------------------------------------------------------------------------------------
def main():
    """ Preprocess sample folder and create config file for snakemake"""

    __doc__ = "Preprocess of the table and create config file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--input_table",
        dest="input_table",
        help="input table containing the sample information",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--input_dict",
        dest="input_dict",
        help="input dictionary containing the feature name \
              conversion from labkey to snakemake allowed names",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--genomes_path",
        dest="genomes_path",
        help="path containing the fasta and gtf files for all organisms",
        required=True)

    parser.add_argument(
        "--multimappers",
        dest="multimappers",
        help="number of mulitmappers allowed",
        required=False,
        type=int,
        metavar='value',
        default=1)

    parser.add_argument(
        "--soft_clip",
        dest="soft_clip",
        help="soft-clipping option of STAR",
        required=False,
        choices=['EndToEnd','Local'],
        default='EndToEnd')

    parser.add_argument(
        "--pass_mode",
        dest="pass_mode",
        help="STAR option pass_mode",
        required=False,
        choices=['None','Basic'],
        default='None')

    parser.add_argument(
        "--libtype",
        dest="libtype",
        help="Library type for salmon",
        required=False,
        default='A')

    parser.add_argument(
        "--config_file",
        dest="config_file",
        help="Configuration file to be used by Snakemake",
        required=False)

    parser.add_argument(
        "--samples_table",
        dest="samples_table",
        help="Table with samples",
        required=True)


    # __________________________________________________________________________________________________________________
    # ------------------------------------------------------------------------------------------------------------------
    # get the arguments
    # ------------------------------------------------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    sys.stdout.write('Reading input file...\n')
    input_table = pd.read_csv(
        options.input_table,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')

    input_dict = pd.read_csv(
        options.input_dict,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    input_dict.set_index('snakemake', inplace=True, drop=True)
    sys.stdout.write('Create snakemake table...\n')
    snakemake_table = pd.DataFrame()
    for index, row in input_table.iterrows():
        if row[input_dict.loc['seqmode', 'labkey']] == 'PAIRED':
            snakemake_table.loc[index, 'seqmode'] = 'paired_end'
        elif row[input_dict.loc['seqmode', 'labkey']] == 'SINGLE':
            snakemake_table.loc[index, 'seqmode'] = 'single_end'

        fq1 = os.path.join(
            row[input_dict.loc['fastq_path', 'labkey']],
            row[input_dict.loc['fq1', 'labkey']])
        snakemake_table.loc[index, 'fq1'] = fq1

        with gzip.open(fq1, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_length = len(record.seq)
                break
        snakemake_table.loc[index, 'index_size'] = read_length
        if read_length <= 50:
            snakemake_table.loc[index, 'kmer'] = 21
        elif read_length > 50:
            snakemake_table.loc[index, 'kmer'] = 31


        snakemake_table.loc[index, 'fq2'] = os.path.join(
            row[input_dict.loc['fastq_path', 'labkey']],
            row[input_dict.loc['fq2', 'labkey']])

        snakemake_table.loc[index, 'fq1_3p'] = row[input_dict.loc['fq1_3p', 'labkey']]
        snakemake_table.loc[index, 'fq1_5p'] = row[input_dict.loc['fq1_5p', 'labkey']]
        snakemake_table.loc[index, 'fq2_3p'] = row[input_dict.loc['fq2_3p', 'labkey']]
        snakemake_table.loc[index, 'fq2_5p'] = row[input_dict.loc['fq2_5p', 'labkey']]

        organism = row[input_dict.loc['organism', 'labkey']].replace(' ', '_').lower()
        snakemake_table.loc[index, 'organism'] = organism
        snakemake_table.loc[index, 'gtf'] = os.path.join(
            options.genomes_path,
            organism,
            'annotation.gtf')
        snakemake_table.loc[index, 'gtf_filtered'] = os.path.join(
            options.genomes_path,
            organism,
            'annotation.gtf')
        snakemake_table.loc[index, 'genome'] = os.path.join(
            options.genomes_path,
            organism,
            'genome.fa')
        snakemake_table.loc[index, 'tr_fasta_filtered'] = os.path.join(
            options.genomes_path,
            organism,
            'transcriptome.fa')

        snakemake_table.loc[index, 'sd'] = row[input_dict.loc['sd', 'labkey']]
        snakemake_table.loc[index, 'mean'] = row[input_dict.loc['mean', 'labkey']]
        snakemake_table.loc[index, 'multimappers'] = options.multimappers
        snakemake_table.loc[index, 'soft_clip'] = options.soft_clip
        snakemake_table.loc[index, 'pass_mode'] = options.pass_mode
        snakemake_table.loc[index, 'libtype'] = options.libtype

        if row[input_dict.loc['mate1_direction', 'labkey']] == 'SENSE':
            snakemake_table.loc[index, 'kallisto_directionality'] = '--fr-stranded'
        elif row[input_dict.loc['mate1_direction', 'labkey']] == 'ANTISENSE':
            snakemake_table.loc[index, 'kallisto_directionality'] = '--rf-stranded'
        else:
            snakemake_table.loc[index, 'kallisto_directionality'] = ''

        if row[input_dict.loc['mate1_direction', 'labkey']] == 'SENSE':
            snakemake_table.loc[index, 'fq1_polya'] = 'AAAAAAAAAAAAAAAAA'
        elif row[input_dict.loc['mate1_direction', 'labkey']] == 'ANTISENSE':
            snakemake_table.loc[index, 'fq1_polya'] = 'TTTTTTTTTTTTTTTTT'
        elif row[input_dict.loc['mate1_direction', 'labkey']] == 'RANDOM':
            snakemake_table.loc[index, 'fq1_polya'] = 'AAAAAAAAAAAAAAAAA'
        else:
            pass

        if row[input_dict.loc['mate2_direction', 'labkey']] == 'SENSE':
            snakemake_table.loc[index, 'fq2_polya'] = 'AAAAAAAAAAAAAAAAA'
        elif row[input_dict.loc['mate2_direction', 'labkey']] == 'ANTISENSE':
            snakemake_table.loc[index, 'fq2_polya'] = 'TTTTTTTTTTTTTTTTT'
        elif row[input_dict.loc['mate2_direction', 'labkey']] == 'RANDOM':
            snakemake_table.loc[index, 'fq2_polya'] = 'AAAAAAAAAAAAAAAAA'
        else:
            pass

    snakemake_table.to_csv(
        options.samples_table,
        sep='\t',
        header=True,
        index=False)

    # Read file and infer read size for sjdbovwerhang
    with open(options.config_file, 'w') as config_file:
        config_file.write('''---
  output_dir: "results"
  local_log: "local_log"
  star_indexes: "star_indexes"
  kallisto_indexes: "kallisto_indexes"
...''')


    sys.stdout.write('Create snakemake table finished successfully...\n')
    sys.stdout.write('Create config file...\n')
    sys.stdout.write('Create config file finished successfully...\n')
    return


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)


