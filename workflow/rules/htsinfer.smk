"""Independent Snakefile for inferring sample specific parameters with HTSinfer."""

import pandas as pd

# set config defaults if not given
try:
    config["records"]
except KeyError:
    config["records"] = 100000
try:
    config["outdir"] 
except KeyError:
    config["outdir"] = os.path.join("tests","input_files")
try:
    config["samples"]
except KeyError:
    config["samples"] = os.path.join("tests","input_files", "samples_in.tsv")
try: 
    config["samples_out"]
except KeyError:
    config["samples_out"] = "samples_htsinfer.tsv"

# global variables
samples = pd.read_csv(config["samples"], header=0, index_col=0, sep = "\t")
OUT_DIR = config["outdir"]
# Write inferred params into new sample table.
SAMPLES_OUT = os.path.join(OUT_DIR,config["samples_out"])



localrules: all, run_htsinfer, htsinfer_to_tsv

rule all:
    input:
        SAMPLES_OUT


# rule run_htsinfer:
#     ''' Infer sample specific parameters'''
#     input:
#         fq1_path = lambda wildcards:
#                 samples.loc[wildcards.sample,"fq1"],
#         fq2_path = lambda wildcards:
#                 samples.loc[wildcards.sample,"fq2"]
#     output:
#         htsinfer_json = os.path.join(OUT_DIR, "htsinfer_{sample}.json")
#     params:
#         records = config["records"]
#         outdir = config["outdir"]
#     threads: 4
#     shell:
#         """
#         htsinfer \
#             --records={params.records} \
#             --output-directory={params.outdir} \
#             --temporary-directory={resources.tmpdir} \
#             --cleanup-regime=KEEP_ALL \
#             --threads={threads} \
#             {input.fq1_path} {input.fq2_path} \
#             > {output.htsinfer_json}
#         """

rule htsinfer_to_tsv:
    '''Write inferred params for all samples to samples.tsv'''
    input:
        jlist = expand(os.path.join(OUT_DIR, 
            "htsinfer_{sample}.json"), 
            sample = samples.index.values),
        samples_in = config["samples"],
        script = os.path.join("workflow","scripts","htsinfer_to_tsv.py")
    output:
        SAMPLES_OUT
    shell:
        '''
        python {input.script} \
            -f {input.jlist} \
            -s {input.samples_in} \
            -o {output}
        '''
        
onsuccess:
    print("Workflow finished, no error.")
onerror:
    print("Ooops... something went wrong")
