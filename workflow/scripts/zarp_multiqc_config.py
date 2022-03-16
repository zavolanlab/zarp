#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Maria Katsantoni, Maciek Bak
# Company: Mihaela Zavolan, Biozentrum, Basel
# This script is part of the Zavolan lab ZARP pipeline.
# In this script the config file used by multiqc
# (https://multiqc.info) is created.
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os


def main():
    """ Create config file for multiqc"""

    __doc__ = "Create config file for multiqc"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--config",
                        help="Output file destination for config",
                        required=True,
                        metavar="FILE",)

    parser.add_argument("--intro-text",
                        dest="intro_text",
                        help="short description at the top of report",
                        metavar="STR")

    parser.add_argument("--custom-logo",
                        dest="custom_logo",
                        default='None',
                        help="Logo path",
                        metavar="FILE")

    parser.add_argument("--url",
                        help="Url of the lab",
                        metavar="STR")

    parser.add_argument("--author-name",
                        dest="author_name",
                        default='None',
                        help="Name of person running this analysis",
                        metavar="STR")

    parser.add_argument("--author-email",
                        dest="author_email",
                        default='None',
                        help="email of person running this analysis",
                        metavar="STR")

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    title = "ZARP"
    subtitle = "RNA-Seq processing pipeline developed by Zavolan Lab"
    logo_title = "ZARP"
    project_type = "Snakemake workflow"
    analysis_type = "RNA-seq"

    intro_text = options.intro_text
    custom_logo = options.custom_logo
    url = options.url
    author_name = options.author_name
    author_email = options.author_email

    config_string = f"""---

title: "{title}"
subtitle: "{subtitle}"
intro_text: "{intro_text}"
custom_logo: "{custom_logo}"
custom_logo_url: "{url}"
custom_logo_title: "{logo_title}"

report_header_info:
  - Project Type: "{project_type}"
  - Analysis Type: "{analysis_type}"
  - Analysis Author: "{author_name}"
  - Contact E-mail: "{author_email}"

module_order:

  - fastqc:
      name: "FastQC: raw reads"
      anchor: "fastqc_raw"
      path_filters:
      - "*/*/fastqc/*/*"

  - cutadapt:
      name: "Cutadapt: adapter removal"
      path_filters:
      - "*/*/*remove_adapters_cutadapt*.stdout.log"

  - cutadapt:
      name: "Cutadapt: polyA tails removal"
      path_filters:
      - "*/*/*remove_polya_cutadapt*.stdout.log"

  - fastqc:
      name: "FastQC: trimmed reads"
      anchor: "fastqc_trimmed"
      path_filters:
      - "*/*/fastqc_trimmed/*/*"

  - star:
      path_filters:
      - "*/*/map_genome/*"

  - ALFA:
      path_filters:
      - "*/*/ALFA/*/*ALFA_feature_counts.tsv"

  - tin-score:
      path_filters:
      - "*/*/TIN/TIN_score.tsv"

  - salmon:
      path_filters:
      - "*/*/*.salmon.*/*"

  - kallisto:
      path_filters:
      - "*/*/*genome_quantification_kallisto*.stderr.log"

  - zpca:
      name: "zpca: salmon | gene expression"
      path_filters:
      - "*/zpca/pca_salmon_genes/*"

  - zpca:
      name: "zpca: salmon | transcript expression"
      path_filters:
      - "*/zpca/pca_salmon_transcripts/*"

  - zpca:
      name: "zpca: kallisto | gene expression"
      path_filters:
      - "*/zpca/pca_kallisto_genes/*"

  - zpca:
      name: "zpca: kallisto | transcript expression"
      path_filters:
      - "*/zpca/pca_kallisto_transcripts/*"

fn_clean_exts:
  - '.fq1'
  - '.gz'
  - '.stdout'
  - '.log'
  - '.stderr'
  - '.fastq'
  - '.bam'
  - '.bai'
  - '.pe'
  - '.se'
  - '.pseudo'
  - '.salmon'
  - '.sam'
..."""

    with open(options.config, "w") as config:
        config.write(config_string)

    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)
