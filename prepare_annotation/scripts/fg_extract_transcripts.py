# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
    import HTSeq
except(Exception):
    raise("[ERROR] HTSeq was not imported properly. Exiting.")
    sys.exit(-1)

try:
    from argparse import ArgumentParser, RawTextHelpFormatter, FileType
except(Exception):
    raise("[ERROR] argparse was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import sys
except(Exception):
    raise("[ERROR] sys was not imported properly. Exiting.")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():

    """ Main function """

    __doc__ = "Filter gtf file based on the transcript biotype."

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)


    parser.add_argument("--gtf",
                      dest = "gtf",
                      help = "Annotation file in GTF format with transcript support level information",
                      required = True,
                      metavar = "FILE")

    parser.add_argument("--out",
                      dest = "out",
                      help = "GTF output file",
                      required = True,
                      metavar="FILE")

    parser.add_argument("--transcript_biotype",
                      dest = "transcript_biotype",
                      help = "Transcript biotype or transcript type to include",
                      required = True,
                      default = "protein_coding,lincRNA")

    parser.add_argument("-v",
                      "--verbose",
                      action = "store_true",
                      dest = "verbose",
                      default = False,
                      required = False,
                      help = "Verbose")

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
        sys.stdout.write("Parsing gtf file and keep %s transcripts\n" % (str(options.transcript_biotype)))

    list_of_transcript_types = str(options.transcript_biotype).split(",")

    # parse gtf file
    gtf = HTSeq.GFF_Reader(options.gtf)
    # open output file
    w = open(options.out, 'w')
    for gtf_line in gtf:
        if "transcript_biotype" in gtf_line.attr:
            if gtf_line.attr['transcript_biotype'] in list_of_transcript_types:
                w.write(gtf_line.get_gff_line())
        if "transcript_type" in gtf_line.attr:
            if gtf_line.attr['transcript_type'] in list_of_transcript_types:
                w.write(gtf_line.get_gff_line())
    w.close()


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
