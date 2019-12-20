#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    __doc__ = "Create heatmap and clustermap based on a table (rows are genes/transcripts, columns are samples)."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--tpm",
        dest="tpm",
        help="TPM table (tsv)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output directory",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose"
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.verbose:
        sys.stdout.write(f"Creating output directory: {options.out} {os.linesep}")

    if not os.path.exists(options.out):
        os.makedirs(options.out)

    if options.verbose:
        sys.stdout.write(f"Reading: {options.tpm} {os.linesep}")
    df = pd.read_csv(options.tpm, header=0, sep="\t")
    df.set_index(["Name"], inplace=True)
    df = df + 1/32
    df = np.log2(df)

    if options.verbose:
        sys.stdout.write(f"Generating heatmap {os.linesep}")

    plt.rcParams["figure.figsize"] = (20,20)
    g = sns.heatmap(df.corr(method="pearson"))
    plt.savefig(os.path.join(options.out, "heatmap.pdf"))

    if options.verbose:
        sys.stdout.write(f"Generating clustermap {os.linesep}")

    plt.rcParams["figure.figsize"] = (20,20)
    g = sns.clustermap(df.corr(method="pearson"))
    plt.savefig(os.path.join(options.out, "clustermap.pdf"))

    if options.verbose:
        sys.stdout.write(f"Done {os.linesep}")

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
