from encodings import utf_8
import os
import json
import logging
import datetime
from jsonschema import ValidationError
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter

# TO DO ############################
# -polyA ? we only need yes or no for 3p and 5p, but where is this handled?
# - index size? We could check the first n reads and take the max. Where to handle this?
# kmer: default 31, where to set?
# organism: has to be translated to the same name as in the annotation files!
# 

def parse_arguments():
    '''
    Parser of command-line arguments.
    '''
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f", 
        "--file_list", 
        help="list of htsinfer output json filepaths.", required=True, nargs="+"
    )
    parser.add_argument(
        "-s", 
        "--samples_in",
        help="path to input samples table", 
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        help="path to output samples table",
        required=True
    )
    return parser


def main():
    # input parameters
    file_list = options.file_list
    samples_df = pd.read_csv(options.samples_in, header=0, index_col=0, sep = "\t")
    outfile = options.output

    params_df = pd.DataFrame()

    # for each sample, load htsinfer json and get params
    for sample in samples_df.index:

        # Get path to htsinfer json for current sample
        filename = [n for n in file_list if sample in n][0]

        logging.debug(f"Filename: {filename}")

        with open(filename, mode='r', encoding="utf-8") as f:
                jparams = json.load(f)

        # Call function to convert the json values into a pd.Series
        tparams = htsinfer_to_zarp(sample,jparams)
        
        # Add to the DataFrame for all samples
        params_df = params_df.append(tparams)
        logging.debug(f"df: {params_df}")


    # Add the new params to samples table:
    # If value is already present, do nothing (prioritize user spec)
    # If value or col missing (missing values must be NaN), fill in from htsinfer.
    
    # replace empty strings
    samples_df = samples_df.replace(r'^\s*$', np.nan, regex=True)
    # replace None
    samples_df =  samples_df.fillna(value=np.nan)
    # combine samples table and inferred params
    samples_df = samples_df.combine_first(params_df)

    # And write new samples table to file
    with open(outfile, mode='w', encoding='utf-8') as o:
        samples_df.to_csv(o,sep="\t",header=True)



def htsinfer_to_zarp(sample,jparams):
    '''Translate htsinfer json output to zarp compatible row.'''  
    tparams = pd.Series(name=sample, dtype="object")

    # seqmode
    f1 = jparams["library_type"]["file_1"]
    f2 = jparams["library_type"]["file_2"]
    rel = jparams["library_type"]["relationship"]
    if  (f1 == "first_mate") and not f2:
        tparams["seqmode"] = "se"
    elif rel == "split_mates":     
        tparams["seqmode"] = "pe"
    else:
        tparams["seqmode"] = None
        logging.warning(f"The sequencing mode could not be determined, as file1 appears to be {f1}, file2 {f2}, and their relationship {rel}.")

    # fq1_3p (same for se and pe)
    f1_3p = jparams["read_layout"]["file_1"]["adapt_3"]
    if f1_3p:
        tparams["fq1_3p"] = f1_3p
    else:
        tparams["fq1_3p"] = "X" * 15
        logging.info("No 3p adapter for fq1 identified.")
            
    # fq2_3p
    f2_3p = jparams["read_layout"]["file_2"]["adapt_3"]
    if f2_3p:
        tparams["fq2_3p"] = f2_3p
    else:
        tparams["fq2_3p"] = "X" * 15

        # info only necessary for pe mode
        if tparams["seqmode"] == "pe":
            logging.info("No 3p adapter for fq2 identified, no adapters will be removed.")

    # organism
    org1 = jparams["library_source"]["file_1"]["short_name"]
    org2 = jparams["library_source"]["file_2"]["short_name"]
    if org1:
        if org2 and (org1 != org2):  
            tparams["organism"] = None
            logging.warning(f"The library source could not be determined, as file1 seems to be derived from {org1}, while file 2 seems to be derived from {org2}.")
        else:
            tparams["organism"] = org1
    elif not org2:
        tparams["organism"] = None
        logging.warning(f"The library source could not be determined, please specify an organism yourself")
    else:
        tparams["organism"] = org2

    # libtype
    f1_o = jparams["read_orientation"]["file_1"]
    rel_o = jparams["read_orientation"]["relationship"]
    if rel_o:
        tparams["libtype"] = rel_o
    elif f1_o:
        tparams["libtype"] = f1_o
    else:
        tparams["libtype"] = None
        logging.warning("The read orientation could not be determined; No values found.")

    # index size
    read_lengths = []
    for i in ["file_1", "file_2"]:
        read_lengths.extend(jparams["library_stats"][i]["read_length"].values())
    
    if read_lengths:
        tparams["index_size"] = int(max([i for i in read_lengths if i]))
    else:
        logging.warning("Read lengths (=index_size) could not be determined")

    # return the pd.Series containing the inferred columns for the current sample
    logging.debug(f"Params: {tparams}")
    return tparams



if __name__ == '__main__':
    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
        # execute the body of the script
        logging.info("Starting script")
        main()
        logging.info("Finished script successfully.")
    
    except Exception as e:
        logging.exception(str(e))
        raise e