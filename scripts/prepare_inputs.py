#!/usr/bin/env python3

"""Create input table and config for Rhea."""

import argparse
from functools import partial
import gzip
import logging
import math
import os
import sys
from typing import Tuple

from Bio import SeqIO
import labkey
import pandas as pd

logger = logging.getLogger(__name__)


def parse_cli_args() -> argparse.Namespace:
    """
    Parses command line arguments.

    :returns: parsed CLI arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "table",
        type=str,
        default=None,
        help="either local file path of input table *or* name of table on "
             "LabKey instance (see 'LabKey API' options below)",
        metavar="TABLE",
    )

    api = parser.add_argument_group("LabKey API")
    api.add_argument(
        "--labkey-domain",
        type=str,
        default=None,
        help="domain of LabKey instance to query; required for obtaining "
             "input table via LabKey API",
        metavar="STR",
    )
    api.add_argument(
        "--labkey-path",
        type=str,
        default=None,
        help="path to LabKey container that includes specified input table; "
             "required for obtaining input table via LabKey API",
        metavar="STR",
    )

    io = parser.add_argument_group("input/output")
    io.add_argument(
        "--input-to-output-mapping",
        type=argparse.FileType('r'),
        default=os.path.join(
            os.path.dirname(__file__),
            'prepare_inputs.dict.tsv',
        ),
        help="lookup table with mappings from input (LabKey or LabKey-like) "
             "to output (Snakemake) table; default: '%(default)s'",
        metavar="FILE",
    )
    io.add_argument(
        "--resources-dir",
        type=str,
        default=os.getcwd(),
        help="path containing the genome resources for all organisms "
             "(default: %(default)s)",
        metavar="DIR",
    )
    io.add_argument(
        "--output-table",
        type=argparse.FileType('w'),
        default="samples.tsv",
        help="output sample table for use in Rhea (default: %(default)s)",
        metavar="FILE",
    )
    io.add_argument(
        "--config-file",
        type=argparse.FileType('w'),
        default="config.yaml",
        help="output Snakemake configuration file for use in Rhea (default: "
             "%(default)s)",
        metavar="FILE",
    )
    io.add_argument(
        "--output-dir",
        type=str,
        default=os.getcwd(),
        help="directory to which Rhea results and logs are to be written "
             "(default: %(default)s)",
        metavar="DIR",
    )
    parser.add_argument(
        "--no-process-paths",
        action="store_true",
        default=False,
        help="do not attempt to create absolute paths in output files",
    )

    behavior = parser.add_argument_group("workflow behavior")
    behavior.add_argument(
        "--trim-polya",
        type=int,
        choices=[True, False],
        default=True,
        help="cutadapt: trim poly(A) tails option (default: %(default)s)",
    )
    behavior.add_argument(
        "--multimappers",
        type=int,
        default=100,
        help="STAR: number of multimappers to report (default: %(default)s)",
        metavar='INT',
    )
    behavior.add_argument(
        "--soft-clip",
        type=str,
        default="EndToEnd",
        help="STAR: soft-clipping option (default: %(default)s)",
        choices=['EndToEnd', 'Local'],
    )
    behavior.add_argument(
        "--pass-mode",
        type=str,
        default="None",
        help="STAR: 2-pass mode option (default: %(default)s)",
        choices=["None", "Basic"],
    )
    behavior.add_argument(
        "--libtype",
        type=str,
        default="A",
        help="Salmon: library type (default: %(default)s)",
        metavar="STR",
    )

    report = parser.add_argument_group("report")
    report.add_argument(
        "--description",
        type=str,
        default="N/A",
        help="short description to be added to the report (default: "
             "%(default)s)",
        metavar="STR",
    )
    report.add_argument(
        "--logo",
        type=argparse.FileType('r'),
        default=None,
        help="path to image file to be added to the report (default: "
             "%(default)s)",
        metavar="FILE",
    )
    report.add_argument(
        "--url",
        type=str,
        default="N/A",
        help="contact URL to be added to the report (default: %(default)s)",
        metavar="STR",
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="print log messages to STDERR",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="print log and debug messages to STDERR",
    )

    args = parser.parse_args()

    if args.logo:
        args.logo.close()
        args.logo = args.logo.name
    else:
        args.logo = ""

    if (args.labkey_domain and not args.labkey_path) or \
       (args.labkey_path and not args.labkey_domain):
        parser.print_help()
        sys.exit(
            "\n[ERROR] Either none or both of '--labkey-domain' and "
            "'--labkey-path' are required."
        )
    return args


def setup_logging(
    logger: logging.Logger,
    verbose: bool = False,
    debug: bool = False,
) -> None:
    """
    Configure logger.

    :param logger: the `logging.Logger` object to configure
    :param verbose: whether `logging.INFO` messages shall be logged
    :param debug: whether `logging.DEBUG` messages shall be logged

    :returns: None
    :raises ?: TODO
    """
    if debug:
        logger.setLevel(logging.DEBUG)
    elif verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        "[%(asctime)-15s: %(levelname)-8s @ %(funcName)s] %(message)s"
    ))
    logger.addHandler(handler)


def fetch_labkey_table(
    domain: str,
    container_path: str,
    query_name: str,
    context_path: str = "labkey",
    schema_name: str = "lists",
) -> pd.DataFrame:
    """
    Export LabKey table as Pandas data frame.

    :param domain: domain of LabKey instance
    :param container_path: path to LabKey container that includes the table of
        interest
    :param query_name: name of LabKey table to export
    :context_path: required by API; usage unclear TODO
    :schema_name: required by API; usage unclear TODO

    :returns: Pandas data frame
    :raises ?: TODO
    """
    server_context = labkey.utils.create_server_context(
        domain=domain,
        container_path=container_path,
        context_path=context_path,
        use_ssl=True,
    )
    results = labkey.query.select_rows(
        server_context=server_context,
        schema_name=schema_name,
        query_name=query_name,
    )
    input_table = pd.DataFrame(results["rows"])
    return input_table


def get_read_length(file: str) -> int:
    """
    Returns read length of first entry of gzipped FASTQ file.

    :param file: path to gzipped FASTQ file

    :returns: read length
    :raises FileNotFoundError: file does not exist
    :raises IsADirectoryError: file is a directory
    :raises OSError: file is not gzipped
    :raises PermissionError: file cannot be read
    :raises ValueError: not a valid FASTQ file
    """
    with gzip.open(file, "rt") as handle:
        return len(next(SeqIO.parse(handle, "fastq")))


def kmer_from_read_length(
    length: int,
    k_max: int = 31,
    k_min: int = 11,
) -> int:
    """
    Given a read length, returns appropriate kmer parameter size for Salmon
    (https://salmon.readthedocs.io/) or similar k-mer-based quantification
    tools.

    References for implementation:
    https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
    https://groups.google.com/d/msg/sailfish-users/fphjX7OIGzY/bMBwlCaZAgAJ

    :param length: length of read in nucleotides
    :param k_max: maximum allowed k-mer size
    :param k_min: minimum allowed k-mer size

    :returns: k_max for l > 2 * k_max, or else the maximum of k and k_min,
        where k is biggest odd integer that fulfills k < l / 2
    """
    k = k_max
    if length < 2 * k_max + 1:
        # ensure kmer is smaller than half of read length
        k = math.floor((length - 1) / 2)
        # ensure kmer is odd
        if not k % 2:
            k -= 1
    if k < k_min:
        k = k_min
    return k


def get_strand_param_kallisto(directionality: str) -> str:
    """
    Returns appropriate strand info parameter for kallisto
    (https://pachterlab.github.io/kallisto/), given a string indicating the
    "directionality" of a sequencing library.

    :param directionality: direction in which library was sequenced (one of
        "SENSE" and "ANTISENSE")

    :returns: appropriate kallisto option for specified directionality; an
        empty string is returned if the directionality value is empty or not
        recognized
    """
    if directionality == "SENSE":
        option = "--fr"
    elif directionality == "ANTISENSE":
        option = "--rf"
    else:
        option = ""
    return option


def get_strand_param_alfa(directionality: str) -> str:
    """
    Returns appropriate strand info parameter for ALFA
    (https://github.com/biocompibens/ALFA), given a string indicating the
    "directionality" of a sequencing library.

    :param directionality: direction in which library was sequenced (one of
        "SENSE" and "ANTISENSE")

    :returns: appropriate ALFA option for specified directionality; an empty
        string is returned if the directionality value is empty or not
        recognized
    """
    if directionality == 'SENSE':
        option = 'fr-firststrand'
    elif directionality == 'ANTISENSE':
        option = 'fr-secondstrand'
    else:
        option = ''
    return option


def get_strand_names_alfa(directionality: str) -> Tuple[str, str]:
    """
    Returns appropriate strand name suffixes for ALFA
    (https://github.com/biocompibens/ALFA), given a string indicating the
    "directionality" of a sequencing library.

    :param directionality: direction in which library was sequenced (one of
        "SENSE" and "ANTISENSE")

    :returns: tuple of ALFA strand name suffixes for two coverage tracks of a
        paired-end sequencing library
    """
    if directionality == "SENSE":
        plus = "str1"
        minus = "str2"
    elif directionality == "ANTISENSE":
        minus = "str1"
        plus = "str2"
    else:
        plus = ""
        minus = ""
    return (plus, minus)


def get_polya_adapter_seqs(directionality: str) -> Tuple[str, str]:
    """
    Returns repeat oligomers for detecting and trimming of poly(A) signals from
    a sequencing library, given a string indicating the library's
    "directionality".

    :param directionality: direction in which library was sequenced (one of
        "SENSE" and "ANTISENSE")

    :returns: tuple of two 15-mers to be used to detect and trim poly(A)
        signals from the 3' and 5' ends of the reads of sequencing library,
        respectively
    """
    if directionality == 'SENSE':
        three = 'AAAAAAAAAAAAAAA'
        five = 'XXXXXXXXXXXXXXX'
    elif directionality == 'ANTISENSE':
        three = 'XXXXXXXXXXXXXXX'
        five = 'TTTTTTTTTTTTTTT'
    else:
        three = 'XXXXXXXXXXXXXXX'
        five = 'XXXXXXXXXXXXXXX'
    return (three, five)


def expand_path(
    *args: str,
    anchor: str = os.getcwd(),
    expand: bool = True,
    no_abs: bool = False,
) -> str:
    """
    Constructs absolute path.

    Not tested with symbolic links.

    :param args: path fragments which will be joined to the anchor from left
        to right
    :param anchor: path relative to which the path fragments in *args shall
        be interpreted; can be absolute or relative; in the latter case, it is
        interpreted relative to the current working directory; if path
        fragments evaluate to absolute path (either before or after expansion),
        the path will be returned without considering the anchor
    :param expand: whether environment variables and user directories (e.g,
        `~`) shall be expanded
    :param join_only: path fragments in args are joined, but no further
        processing is done

    :returns: absolute path
    """
    suffix = os.path.join(*args)
    if no_abs:
        return suffix
    if os.path.isabs(suffix):
        return os.path.normpath(suffix)
    if expand:
        suffix = os.path.expanduser(
            os.path.expandvars(
                suffix
            )
        )
    if os.path.isabs(suffix):
        return os.path.normpath(suffix)
    anchor = os.path.expanduser(
        os.path.expandvars(
            anchor
        )
    )
    path = os.path.join(anchor, suffix)
    return os.path.normpath(path)


def main(args):
    """
    Create input table and config for Rhea.
    """
    setup_logging(
        logger=logger,
        verbose=args.verbose,
        debug=args.debug,
    )

    # get input table from LabKey or CLI
    if args.labkey_domain:
        logger.info(
            f"Fetching input table from LabKey instance "
            "'{args.labkey_domain}'..."
        )
        input_table = fetch_labkey_table(
            domain=args.labkey_domain,
            container_path=args.labkey_path,
            query_name=args.table,
        )
        labkey_table = expand_path(
            '.'.join([args.output_table.name, "labkey"])
        )
        input_table.to_csv(
            labkey_table,
            sep='\t',
            index=False,
        )
        from_api = True
    else:
        logger.info(f"Reading input table from file '{args.table}'...")
        input_table = pd.read_csv(
            args.table,
            header=0,
            sep='\t',
            index_col=None,
            comment='#',
            engine='python',
        )
        from_api = False

    # get LabKey to Snakemake sample table field mappings
    input_dict = pd.read_csv(
        args.input_to_output_mapping,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python',
    )
    args.input_to_output_mapping.close()
    input_dict.set_index('snakemake', inplace=True, drop=True)

    # create Snakemake table
    logger.info("Creating Snakemake input table...")
    snakemake_table = pd.DataFrame()

    for index, row in input_table.iterrows():

        # extract data from LabKey-like table
        lk_replicate_name = row[input_dict.loc['replicate_name', 'labkey']]
        lk_condition = row[input_dict.loc['condition', 'labkey']]
        lk_seqmode = row[input_dict.loc['seqmode', 'labkey']]
        lk_fastq_path = row[input_dict.loc['fastq_path', 'labkey']]
        lk_fq1 = row[input_dict.loc['fq1', 'labkey']]
        lk_fq2 = row[input_dict.loc['fq2', 'labkey']]
        lk_fq1_3p = row[input_dict.loc['fq1_3p', 'labkey']]
        lk_fq1_5p = row[input_dict.loc['fq1_5p', 'labkey']]
        lk_fq2_3p = row[input_dict.loc['fq2_3p', 'labkey']]
        lk_fq2_5p = row[input_dict.loc['fq2_5p', 'labkey']]
        lk_organism = row[input_dict.loc['organism', 'labkey']]
        lk_sd = row[input_dict.loc['sd', 'labkey']]
        lk_mean = row[input_dict.loc['mean', 'labkey']]
        lk_mate1_direction = row[input_dict.loc['mate1_direction', 'labkey']]
        lk_mate2_direction = row[input_dict.loc['mate2_direction', 'labkey']]

        # extract, infer or convert to Snakemake input format
        if from_api and not os.path.isabs(lk_fastq_path):
            anchor = os.getcwd()
            logger.warning(
                f"[WARNING] Don't know how to interpret relative paths "
                "inside LabKey table. Trying with current working directory "
                f"'{anchor}' as an anchor, but it may be better to use"
                "absolute paths wherever possible..."
            )
        else:
            anchor = os.path.abspath(os.path.dirname(args.table))
        sample = "_".join([lk_replicate_name, lk_condition])
        if lk_seqmode == 'PAIRED':
            seqmode = 'pe'
            fq2 = expand_path(
                lk_fastq_path,
                lk_fq2,
                anchor=anchor,
            )
        elif lk_seqmode == 'SINGLE':
            seqmode = 'se'
            fq2 = "XXXXXXXXXXXXXXX"
        else:
            logger.error(
                f"[ERROR] Illegal sequencing mode '{lk_seqmode}' in row "
                f"{index+1}."
            )
            sys.exit("Execution aborted.")
        fq1 = expand_path(
            lk_fastq_path,
            lk_fq1,
            anchor=anchor,
        )
        read_length = get_read_length(fq1)
        index_size = read_length - 1
        kmer = kmer_from_read_length(read_length)
        fq1_3p = lk_fq1_3p
        fq1_5p = lk_fq1_5p
        fq2_3p = lk_fq2_3p
        fq2_5p = lk_fq2_5p
        organism = lk_organism.replace(' ', '_').lower()
        gtf = expand_path(
            args.resources_dir,
            organism,
            'annotation.gtf',
        )
        genome = expand_path(
            args.resources_dir,
            organism,
            'genome.fa',
        )
        sd = lk_sd
        mean = lk_mean
        fq1_polya_3p, fq1_polya_5p = get_polya_adapter_seqs(lk_mate1_direction)
        fq2_polya_3p, fq2_polya_5p = get_polya_adapter_seqs(lk_mate2_direction)
        kallisto_directionality = get_strand_param_kallisto(lk_mate1_direction)
        alfa_directionality = get_strand_param_alfa(lk_mate1_direction)
        alfa_plus, alfa_minus = get_strand_names_alfa(lk_mate1_direction)

        # construct row in Snakemake input table
        snakemake_table.loc[index, 'sample'] = sample
        snakemake_table.loc[index, 'seqmode'] = seqmode
        snakemake_table.loc[index, 'fq1'] = fq1
        snakemake_table.loc[index, 'fq2'] = fq2
        snakemake_table.loc[index, 'index_size'] = index_size
        snakemake_table.loc[index, 'kmer'] = kmer
        snakemake_table.loc[index, 'fq1_3p'] = fq1_3p
        snakemake_table.loc[index, 'fq1_5p'] = fq1_5p
        snakemake_table.loc[index, 'fq2_3p'] = fq2_3p
        snakemake_table.loc[index, 'fq2_5p'] = fq2_5p
        snakemake_table.loc[index, 'organism'] = organism
        snakemake_table.loc[index, 'gtf'] = gtf
        snakemake_table.loc[index, 'genome'] = genome
        snakemake_table.loc[index, 'sd'] = sd
        snakemake_table.loc[index, 'mean'] = mean
        snakemake_table.loc[index, 'kallisto_directionality'] = \
            kallisto_directionality
        snakemake_table.loc[index, 'alfa_directionality'] = alfa_directionality
        snakemake_table.loc[index, 'alfa_plus'] = alfa_plus
        snakemake_table.loc[index, 'alfa_minus'] = alfa_minus

        # add CLI argument-dependent parameters
        snakemake_table.loc[index, 'multimappers'] = args.multimappers
        snakemake_table.loc[index, 'soft_clip'] = args.soft_clip
        snakemake_table.loc[index, 'pass_mode'] = args.pass_mode
        snakemake_table.loc[index, 'libtype'] = args.libtype
        if args.trim_polya is True:
            snakemake_table.loc[index, 'fq1_polya_3p'] = fq1_polya_3p
            snakemake_table.loc[index, 'fq1_polya_5p'] = fq1_polya_5p
            snakemake_table.loc[index, 'fq2_polya_3p'] = fq2_polya_3p
            snakemake_table.loc[index, 'fq2_polya_5p'] = fq2_polya_5p

    # adjust sample table format
    snakemake_table.fillna('XXXXXXXXXXXXXXX', inplace=True)
    snakemake_table = snakemake_table.astype(
        {
            "sd": int,
            "mean": int,
            "multimappers": int,
            "kmer": int,
            "index_size": int,
        }
    )

    # write Snakemake sample table
    logger.info("Writing Snakemake input table...")
    snakemake_table.to_csv(
        args.output_table,
        sep='\t',
        header=True,
        index=False)
    args.output_table.close()

    # compile entries for Snakemake config file
    logger.info("Creating Snakemake config file...")
    results_dir = expand_path(
        args.output_dir,
        "results",
    )
    log_dir = expand_path(
        args.output_dir,
        "logs",
    )
    kallisto_indexes = expand_path(
        results_dir,
        "kallisto_indexes",
    )
    salmon_indexes = expand_path(
        results_dir,
        "salmon_indexes",
    )
    star_indexes = expand_path(
        results_dir,
        "star_indexes",
    )
    alfa_indexes = expand_path(
        results_dir,
        "alfa_indexes",
    )

    # write Snakemake config file
    logger.info("Writing Snakemake config file...")
    config_file_content = f'''---
  samples: "{expand_path(args.output_table.name)}"
  output_dir: "{results_dir}"
  log_dir: "{log_dir}"
  kallisto_indexes: "{kallisto_indexes}"
  salmon_indexes: "{salmon_indexes}"
  star_indexes: "{star_indexes}"
  alfa_indexes: "{alfa_indexes}"
  report_description: "{args.description}"
  report_logo: "{args.logo}"
  report_url: "{args.url}"
...
'''
    args.config_file.write(config_file_content)
    args.config_file.close()


if __name__ == '__main__':
    args = parse_cli_args()

    # Set default according to CLI arg
    expand_path = partial(expand_path, no_abs=args.no_process_paths)  # type: ignore

    main(args)
    logger.info("Program completed successfully.")
    sys.exit(0)
