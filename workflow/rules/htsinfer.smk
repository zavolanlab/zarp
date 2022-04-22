"""Independent Snakefile for inferring sample specific parameters with HTSinfer."""

import os
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
LOG_DIR = os.path.join(OUT_DIR,"logs")
CLUSTER_LOG = os.path.join(LOG_DIR,"cluster_logs")
# Write inferred params into new sample table.
SAMPLES_OUT = os.path.join(OUT_DIR,config["samples_out"])



localrules: all, htsinfer_to_tsv

rule all:
    input:
        SAMPLES_OUT

current_rule = "run_htsinfer"
rule run_htsinfer:
    ''' Placeholder to test wiring.
    for real rule, remove input.mock and replace shell call with
    
    htsinfer \
        --records={params.records} \
        --output-directory={params.outdir} \
        --temporary-directory={resources.tmpdir} \
        --cleanup-regime=KEEP_ALL \
        --threads={threads} \
        {input.fq1_path} {params.fq2_path} \
        > {output.htsinfer_json}
    
    
    '''
    input:
        fq1_path = lambda wildcards:
                samples.loc[wildcards.sample, "fq1"],
        mock = os.path.join(OUT_DIR, "{sample}.json")
    output:
        htsinfer_json = os.path.join(OUT_DIR, "htsinfer_{sample}.json")
    params:
        fq2_path = lambda wildcards: 
                (samples.loc[wildcards.sample,"fq2"] if samples.loc[wildcards.sample,"fq2"] else ""),
        records = config["records"],
        outdir = OUT_DIR,
        cluster_log_path = CLUSTER_LOG
    threads: 4
    singularity: 
        "docker://quay.io/biocontainers/htsinfer:0.1"
    conda: 
        os.path.join(workflow.basedir, "envs", "htsinfer.yaml")
    log:
        stderr = os.path.join(
            LOG_DIR,
            "{sample}",
            current_rule + ".stderr.log"),
        stdout = os.path.join(
            LOG_DIR,
            "{sample}",
            current_rule + ".stdout.log")
    shell:
        '''
        cp {input.mock} {output.htsinfer_json}
        echo "fq2"; echo {params.fq2_path}
        '''

current_rule = "htsinfer_to_tsv"
rule htsinfer_to_tsv:
    '''Write inferred params for all samples to samples.tsv'''
    input:
        jlist = expand(os.path.join(OUT_DIR, 
            "htsinfer_{sample}.json"), 
            sample = samples.index.tolist()),
        samples_in = config["samples"],
        script = os.path.join("workflow","scripts","htsinfer_to_tsv.py")
    output:
        SAMPLES_OUT
    threads: 4
    log:
        stderr = os.path.join(
            LOG_DIR,
            current_rule + ".stderr.log"),
        stdout = os.path.join(
            LOG_DIR,
            current_rule + ".stdout.log")
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
