from argparse import ArgumentParser, RawTextHelpFormatter
import json
import logging
import os
import sys

from htsinfer.models import Results
import numpy as np
import pandas as pd

logging.basicConfig(
            level=logging.INFO, 
            format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
LOGGER = logging.getLogger(__name__)



def parse_arguments():
    '''
    Parser of command-line arguments.
    '''
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f", 
        "--file-list", 
        help="list of htsinfer output json filepaths. ATTENTION: the sample name must be part of the filename.", required=True, nargs="+"
    )
    parser.add_argument(
        "-s", 
        "--samples-in",
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
    # replace empty strings
    samples_df = samples_df.replace(r'^\s*$', np.nan, regex=True)
    # replace None
    samples_df =  samples_df.fillna(value=np.nan)

    outfile = options.output

    params_df = pd.DataFrame(columns=[
        "fq1_3p",
        "fq2_3p",
        "index_size",
        "libtype",
        "organism",
        "seqmode"
        ])

    # for each sample, load htsinfer json and get params
    for sample in samples_df.index:

        # Get path to htsinfer json for current sample
        filename = "htsinfer_" + sample + ".json"
        f = [n for n in file_list if filename == os.path.basename(n)][0]

        LOGGER.debug(f"Filename: {f}")

        with open(f, mode='r', encoding="utf-8") as f:
            try:
                jparams = Results(**json.load(f))
            except Exception as exc:
                LOGGER.error(f"Results are not valid htsinfer format.")
                raise ValueError from exc

        # Call function to convert the json values into a pd.Series
        tparams, e_flag, swap_paths = htsinfer_to_zarp(sample,jparams, samples_df)
        
        # Add to the DataFrame for all samples
        params_df.loc[sample] = tparams

        # Need to swap filepaths?
        if swap_paths is True:
            fq2 = samples_df.loc[sample, "fq1"]
            samples_df.loc[sample, "fq1"] = samples_df.loc[sample, "fq2"]
            samples_df.loc[sample, "fq2"] = fq2


    LOGGER.debug(f"params_df: {params_df}")

    # Add the new params to samples table:
    # If value is already present, do nothing (prioritize user spec)
    # If value or col missing (missing values must be NaN), fill in from htsinfer.
    samples_df = samples_df.combine_first(params_df)

    # And write new samples table to file
    with open(outfile, mode='w', encoding='utf-8') as o:
        samples_df.to_csv(o,sep="\t",header=True)

    return e_flag


def should_i_flag(df, sample, param):
    '''Only RAISE error if user hasn't specified value either'''

    try:
        user_param = df.loc[sample,param]
    except KeyError:
        user_param = np.nan
    if user_param is np.nan:
        e_flag = True

    return e_flag


def htsinfer_to_zarp(sample,jparams, samples_df):
    '''Translate htsinfer json output to zarp compatible row.'''  

    # Flag errors
    e_flag = False
    # need to swap filepaths?
    swap_paths = False

    # save inferred params as sample row
    tparams = pd.Series(name=sample, dtype="object")

    # seqmode
    # NOTE: library_type contents have to be checked explicitly, as the corresponding objects in htsinfer model are optional. (so if they are not present, '.value' will cause an error) 
    if jparams.library_type.file_1 is not None:
        f1 = jparams.library_type.file_1.value
    else:
        f1 = None

    if jparams.library_type.file_2 is not None:
        f2 = jparams.library_type.file_2.value
    else:
        f2 = None

    if jparams.library_type.relationship is not None:
        rel = jparams.library_type.relationship.value
    else:
        rel = None
        
    if f1 is not None and f2 is None:
        tparams["seqmode"] = "se"
    elif rel == "split_mates":     
        tparams["seqmode"] = "pe"
        if f1 == "second_mate":
            swap_paths = True
    else:
        tparams["seqmode"] = None
        LOGGER.error(f"The sequencing mode could not be determined, as file1 appears to be {f1}, file2 {f2}, and their relationship {rel}.")
        
        e_flag = should_i_flag(samples_df, sample, "seqmode")
        

    # fq1_3p (same for se and pe)
    f1_3p = jparams.read_layout.file_1.adapt_3
    tparams["fq1_3p"] = "X" * 15 if f1_3p is None else f1_3p    
    if f1_3p is None:
        LOGGER.warning("No 3p adapter for fq1 identified.")
            
    # fq2_3p
    f2_3p = jparams.read_layout.file_2.adapt_3
    tparams["fq2_3p"] = "X" * 15 if f2_3p is None else f2_3p
    if f2_3p is None:
        # info only necessary for pe mode
        if tparams["seqmode"] == "pe":
            LOGGER.warning("No 3p adapter for fq2 identified, no adapters will be removed.")

    # organism
    org1 = jparams.library_source.file_1.short_name
    org2 = jparams.library_source.file_2.short_name
    tparams["organism"] = None

    # source could not be inferred
    if org1 is None and org2 is None:
        LOGGER.error(f"The library source could not be determined, please specify an organism yourself")
        
        e_flag = should_i_flag(samples_df, sample, "organism")

    # source inferred for single-end library or one mate of a paired-end library (used for both)
    elif None in (org1, org2):
        tparams["organism"] = org2 if org1 is None else org1

    # same source for both mates of a paired-end library
    elif org1 == org2:
        tparams["organism"] = org1

    # different sources inferred for both mates
    else:
        LOGGER.error(f"The library source could not be determined, as file1 seems to be derived from {org1}, while file 2 seems to be derived from {org2}.")
        
        e_flag = should_i_flag(samples_df, sample, "organism")

    # libtype
    # NOTE: read_orientation contents have to be checked explicitly, as the corresponding objects in htsinfer model are optional. (so if None, .value will cause an error) 
    if jparams.read_orientation.file_1 is not None:
        f1_o = jparams.read_orientation.file_1.value
    else:
        f1_o = None

    if jparams.read_orientation.relationship is not None:
        rel_o = jparams.read_orientation.relationship.value
    else:
        rel_o = None

    if rel_o is not None:
        tparams["libtype"] = rel_o
    elif f1_o is not None:
        tparams["libtype"] = f1_o
    else:
        tparams["libtype"] = None
        LOGGER.error("The read orientation could not be determined; No values found.")
        
        e_flag = should_i_flag(samples_df, sample, "libtype")

    # index size
    read_lengths = []
    read_lengths.append(jparams.library_stats.file_1.read_length.max)
    read_lengths.append(jparams.library_stats.file_2.read_length.max)
    if (read_lengths is not None) and (len(read_lengths) != 0):
        tparams["index_size"] = max([int(i) for i in read_lengths if i is not None])
    else:
        LOGGER.error("Read lengths (=index_size) could not be determined")
        
        e_flag = should_i_flag(samples_df, sample, "index_size")

    # return the pd.Series containing the inferred columns for the current sample, the error flag, and filepath swap instruction
    LOGGER.debug(f"Params: {tparams}")
    return tparams, e_flag, swap_paths



if __name__ == '__main__':
    # parse the command-line arguments
    options = parse_arguments().parse_args()
    
    # execute the body of the script
    LOGGER.info("Starting script")
    e_flag = main()

    if e_flag is False:
        LOGGER.info("Finished script successfully.")
    else:
        LOGGER.error("Finished script BUT one or more required parameters could not be inferred, please check logs.")
        sys.exit(1)
