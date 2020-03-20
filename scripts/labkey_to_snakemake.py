#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Katsantoni Maria, Christina Herrmann
# Company: Mihaela Zavolan, Biozentrum, Basel
# This script is part of the Zavolan lab quantification pipeline, which is used
# for analysing RNA-seq data. The table is provided by labkey as a csv file.
# If the user provides their own table the table should contain the following
# columns:
# -----------------------------------------------------------------------------

import sys
import gzip
import labkey
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO
from csv import writer
from pathlib import Path
# (avoids long lines in filter definitions)
from labkey.query import QueryFilter


def main():
    """ Preprocess sample folder and create config file for snakemake"""

    __doc__ = "Preprocess of labkey table and create " + \
              "config file and sample table."

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("genomes_path",
                        help="Path containing the FASTA and GTF " +
                        " files for all organisms",
                        metavar="GENOMES PATH")

    parser.add_argument("--input-table",
                        type=str,
                        default=None,
                        help="Input table in LabKey format " +
                        "containing the sample information;" +
                        "\nexactly one '--input-table' and " +
                        "'--remote' is required.",
                        metavar="FILE")

    parser.add_argument("--remote",
                        action="store_true",
                        help="Fetch LabKey table via API; exactly one of " +
                        "'--input-table' and" +
                        "\n'--remote' is required.")

    parser.add_argument("--project-name",
                        help="Name of LabKey project containing table " +
                        " '--table-name'; required" +
                        "\nif '--remote' is specified.",
                        metavar="STR")

    parser.add_argument("--table-name",
                        help="Name of LabKey table; required if '--remote'" +
                        " is specified.",
                        metavar="STR")

    parser.add_argument("--input-dict",
                        help="Input dictionary containing the feature name " +
                        "conversion from LabKey to Snakemake;" +
                        "default: '%(default)s'",
                        default=os.path.join(
                            os.path.dirname(__file__),
                            'labkey_to_snakemake.dict.tsv'),
                        metavar="FILE")

    parser.add_argument("--samples-table",
                        help="Output table compatible to snakemake;" +
                        "default: '%(default)s'",
                        default='samples.tsv',
                        metavar="FILE")

    parser.add_argument("--trim_polya",
                        type=int,
                        choices=[True, False],
                        default=True,
                        help="Trim poly-As option")

    parser.add_argument("--multimappers",
                        type=int,
                        default=100,
                        help="Number of allowed multimappers",
                        metavar='INT')

    parser.add_argument("--soft-clip",
                        choices=['EndToEnd', 'Local'],
                        default='EndToEnd',
                        help="Soft-clipping option for STAR")

    parser.add_argument("--pass-mode",
                        choices=['None', 'Basic'],
                        default='None',
                        help="2-pass mode option for STAR")

    parser.add_argument("--libtype",
                        default='A',
                        help="Library type for salmon",
                        metavar="STR")

    parser.add_argument("--config-file",
                        help="Configuration file to be used by Snakemake")

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.remote and options.input_table:
        parser.print_help()
        print(
            "\n[ERROR] Options '--input-table' and ",
            "'--remote' are mutually exclusive.")
        sys.exit(1)

    if not options.remote and not options.input_table:
        parser.print_help()
        print("\n[ERROR] At least one of '--input-table' ",
              "and '--remote' is required.")
        sys.exit(1)

    if options.remote and not options.project_name:
        parser.print_help()
        print(
            "\n[ERROR] If option '--remote' is specified, ",
            "option '--project-name' is required.")
        sys.exit(1)

    if options.remote and not options.table_name:
        parser.print_help()
        print(
            "\n[ERROR] If option '--remote' is specified, ",
            "option '--table-name' is required.")
        sys.exit(1)

    sys.stdout.write('Reading input file...\n')

    if options.remote is True:
        input_table = api_fetch_labkey_table(
            project_name=options.project_name,
            query_name=options.table_name)
        input_table.to_csv(options.input_table, sep='\t', index=False)
    else:
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
        snakemake_table.loc[index, 'sample'] = row[
            input_dict.loc['replicate_name', 'labkey']] + "_" + row[
            input_dict.loc['condition', 'labkey']]
        if row[input_dict.loc['seqmode', 'labkey']] == 'PAIRED':
            snakemake_table.loc[index, 'seqmode'] = 'paired_end'
        elif row[input_dict.loc['seqmode', 'labkey']] == 'SINGLE':
            snakemake_table.loc[index, 'seqmode'] = 'single_end'

        fq1 = os.path.join(
            row[input_dict.loc['fastq_path', 'labkey']],
            row[input_dict.loc['fq1', 'labkey']])

        snakemake_table.loc[index, 'fq1'] = fq1
        read_length = get_read_length(fq1)
        snakemake_table.loc[index, 'index_size'] = read_length
        snakemake_table.loc[index, 'kmer'] = infer_kmer_length(read_length)
        snakemake_table.loc[index, 'fq1_3p'] = row[
            input_dict.loc['fq1_3p', 'labkey']]
        snakemake_table.loc[index, 'fq1_5p'] = row[
            input_dict.loc['fq1_5p', 'labkey']]

        organism = row[input_dict.loc['organism', 'labkey']].replace(
            ' ', '_').lower()
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

        snakemake_table.loc[index, 'sd'] = row[
            input_dict.loc['sd', 'labkey']]
        snakemake_table.loc[index, 'mean'] = row[
            input_dict.loc['mean', 'labkey']]
        snakemake_table.loc[index, 'multimappers'] = options.multimappers
        snakemake_table.loc[index, 'soft_clip'] = options.soft_clip
        snakemake_table.loc[index, 'pass_mode'] = options.pass_mode
        snakemake_table.loc[index, 'libtype'] = options.libtype

        if options.trim_polya is True:
            fq1_polya_3p, fq1_polya_5p = trim_polya(
                row[input_dict.loc['mate1_direction', 'labkey']])
            snakemake_table.loc[index, 'fq1_polya_3p'] = fq1_polya_3p
            snakemake_table.loc[index, 'fq1_polya_5p'] = fq1_polya_5p

        snakemake_table.loc[index, 'kallisto_directionality'] = \
            get_kallisto_directionality(
                row[input_dict.loc['mate1_direction', 'labkey']])

        if row[input_dict.loc['seqmode', 'labkey']] == 'PAIRED':
            fq2 = os.path.join(
                row[input_dict.loc['fastq_path', 'labkey']],
                row[input_dict.loc['fq2', 'labkey']])
            snakemake_table.loc[index, 'fq2'] = fq2

            snakemake_table.loc[index, 'fq2_3p'] = row[
                input_dict.loc['fq2_3p', 'labkey']]
            snakemake_table.loc[index, 'fq2_5p'] = row[
                input_dict.loc['fq2_5p', 'labkey']]

            if options.trim_polya is True:
                fq2_polya_3p, fq2_polya_5p = trim_polya(
                    row[input_dict.loc['mate2_direction', 'labkey']])
                snakemake_table.loc[index, 'fq2_polya_3p'] = fq2_polya_3p
                snakemake_table.loc[index, 'fq2_polya_5p'] = fq2_polya_5p

    snakemake_table.fillna('XXXXXXXXXXXXX', inplace=True)
    snakemake_table = snakemake_table.astype(
        {
            "sd": int,
            "mean": int,
            "multimappers": int,
            "kmer": int,
            "index_size": int,
        }
    )
    snakemake_table.to_csv(
        options.samples_table,
        sep='\t',
        header=True,
        index=False)

    # Read file and infer read size for sjdbovwerhang
    with open(options.config_file, 'w') as config_file:
        config_file.write('''---
  samples: "''' + options.samples_table + '''"
  output_dir: "results/"
  log_dir: "logs/"
  kallisto_indexes: "results/kallisto_indexes/"
  salmon_indexes: "results/salmon_indexes/"
  star_indexes: "results/star_indexes/"
  alfa_indexes: "results/alfa_indexes/"
...''')

    sys.stdout.write('Create snakemake table finished successfully...\n')
    sys.stdout.write('Create config file...\n')
    sys.stdout.write('Create config file finished successfully...\n')
    return


def api_fetch_labkey_table(project_name=None, query_name=None):
    group_path = os.path.join('/Zavolan Group', project_name)
    server_context = labkey.utils.create_server_context(
        'labkey.scicore.unibas.ch', group_path, 'labkey', use_ssl=True)
    schema_name = "lists"
    results = labkey.query.select_rows(server_context, schema_name, query_name)
    input_table = pd.DataFrame(results["rows"])
    return input_table


def get_read_length(filename):
    with gzip.open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_length = len(record.seq)
            break
    return read_length


def infer_kmer_length(read_length):
    if read_length <= 50:
            kmer = 21
    elif read_length > 50:
        kmer = 31
    return kmer


def get_kallisto_directionality(directionality):
    if directionality == 'SENSE':
        final_direction = '--fr'
    elif directionality == 'ANTISENSE':
        final_direction = '--rf'
    else:
        final_direction = ''
    return final_direction


def trim_polya(sense):
    if sense == 'SENSE':
        polya_3p = 'AAAAAAAAAAAAAAAAA'
        polya_5p = 'XXXXXXXXXXXXXXXXX'
    elif sense == 'ANTISENSE':
        polya_3p = 'XXXXXXXXXXXXXXXXX'
        polya_5p = 'TTTTTTTTTTTTTTTTT'
    else:
        polya_3p = 'XXXXXXXXXXXXXXXXX'
        polya_5p = 'XXXXXXXXXXXXXXXXX'
    return polya_3p, polya_5p


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
