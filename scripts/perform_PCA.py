#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import pandas as pd
import numpy as np
import random as rd
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    __doc__ = "Perform PCA based on a table (rows are genes/transcripts, columns are samples)."

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

    if options.verbose:
        sys.stdout.write(f"Performing PCA {os.linesep}")

    scaled_data = preprocessing.scale(df.T)

    pca = PCA()
    pca.fit(scaled_data)
    pca_data = pca.transform(scaled_data)

    per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
    labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
    plt.figure()
    plt.rcParams["figure.figsize"] = (20,20)
    plt.bar(x=range(1,len(per_var)+1), height=per_var, tick_label=labels)    
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component')
    plt.title('Scree Plot')
    plt.xticks(rotation='vertical')
    plt.savefig(os.path.join(options.out, "scree_plot.pdf"))
    plt.close()
    
    pca_df = pd.DataFrame(pca_data, index=[*df], columns=labels)

    if options.verbose:
        sys.stdout.write(f"Generating plots in: {options.out} {os.linesep}")

    # -------------------------------------------------------------------------
    # PCA 1st and 2nd component
    # -------------------------------------------------------------------------

    plt.figure()
    plt.rcParams["figure.figsize"] = (20,20)
    plt.scatter(pca_df.PC1, pca_df.PC2)
    plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=22)
    plt.ylabel('PC2 - {0}%'.format(per_var[1]), fontsize=22)
    plt.title('PCA components 1 & 2', fontsize=22)
    for sample in pca_df.index:
        plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))
    plt.savefig(os.path.join(options.out, "PCA_1_2.pdf"))
    plt.close()

    # -------------------------------------------------------------------------
    # PCA 2nd and 3rd component
    # -------------------------------------------------------------------------

    plt.figure()
    plt.rcParams["figure.figsize"] = (20,20)
    plt.scatter(pca_df.PC2, pca_df.PC3)
    plt.xlabel('PC2 - {0}%'.format(per_var[1]), fontsize=22)
    plt.ylabel('PC3 - {0}%'.format(per_var[2]), fontsize=22)
    plt.title('PCA components 2 & 3', fontsize=22)
    for sample in pca_df.index:
        plt.annotate(sample, (pca_df.PC2.loc[sample], pca_df.PC3.loc[sample]))
    plt.savefig(os.path.join(options.out, "PCA_2_3.pdf"))
    plt.close()


    # -------------------------------------------------------------------------
    # PCA 1st and 3rd component
    # -------------------------------------------------------------------------

    plt.figure()
    plt.rcParams["figure.figsize"] = (20,20)
    plt.scatter(pca_df.PC1, pca_df.PC3)
    plt.xlabel('PC1 - {0}%'.format(per_var[0]), fontsize=22)
    plt.ylabel('PC3 - {0}%'.format(per_var[2]), fontsize=22)
    plt.title('PCA components 1 & 3', fontsize=22)
    for sample in pca_df.index:
        plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC3.loc[sample]))
    plt.savefig(os.path.join(options.out, "PCA_1_3.pdf"))
    plt.close()
    

    # -------------------------------------------------------------------------
    # PCA 3D
    # -------------------------------------------------------------------------

    labels = []
    pc1 = []
    pc2 = []
    pc3 = []
    for i, row in pca_df.iterrows():
        labels.append(i)
        pc1.append(row["PC1"])
        pc2.append(row["PC2"])
        pc3.append(row["PC3"])

    fig = plt.figure(figsize=(20,20))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pc1, pc2, pc3)

    ax.set_xlabel('PC1 - {0}%'.format(per_var[0]))
    ax.set_ylabel('PC2 - {0}%'.format(per_var[1]))
    ax.set_zlabel('PC3 - {0}%'.format(per_var[2]))
    for label, x, y, z in zip(labels, pc1, pc2, pc3):
        ax.text(x, y, z, label)
    plt.savefig(os.path.join(options.out, "PCA_3D.pdf"))

    # -------------------------------------------------------------------------
    # Write loading scores
    # -------------------------------------------------------------------------

    if options.verbose:
        sys.stdout.write(f"Writing loading scores in: {options.out} {os.linesep}")

    loading_scores = pd.Series(pca.components_[0], index=list(df.index))
    sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
    sorted_loading_scores.to_csv(os.path.join(options.out, "loading_scores.tsv"), header=False, sep="\t")

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
