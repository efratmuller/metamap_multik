# snakemake workflow for mapping metagenomes against 2 reference databases (herein refered to as database A and database B)

import os

# preparing files
configfile: "config.yml"

INPUT_FILE = config["input"]
OUTPUT_DIR = config["output"]
REF_NAME_A = config["database_a"]["db_name"]
REF_NAME_B = config["database_b"]["db_name"]

os.system("chmod -R +x tools")

# parse the list of samples and paths to fastqs
samp2path = {}
with open(INPUT_FILE) as f:
    for line in f:
        cols = line.strip().split()
        fwd = cols[0]
        rev = cols[1]
        sample = os.path.basename(cols[0]).split("_1.fastq")[0]
        samp2path[sample] = [fwd, rev]

# create folders required for later
for sample in samp2path.keys():
    if not os.path.exists(OUTPUT_DIR+"/mapping/"+sample+"/logs"):
        os.makedirs(OUTPUT_DIR+"/mapping/"+sample+"/logs")
if not os.path.exists(OUTPUT_DIR+"/summary/logs"):
    os.makedirs(OUTPUT_DIR+"/summary/logs")

# rule that specifies the final expected output files
rule all:
    input:
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_A+"_cov-est.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_A+"_cov-exp.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_A+"_counts-total.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_A+"_counts-unique.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_B+"_cov-est.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_B+"_cov-exp.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_B+"_counts-total.csv",
        OUTPUT_DIR+"/summary/bwa_"+REF_NAME_B+"_counts-unique.csv"

# map metagenomes to catalog A
rule map2ref_a:
    input:
        fwd = lambda wildcards: samp2path[wildcards.sample][0],
        rev = lambda wildcards: samp2path[wildcards.sample][1]
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_A+"_total.tab",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_A+"_unique.tab"
    params:
        outpref = OUTPUT_DIR+"/mapping/{sample}/{sample}_" + REF_NAME_A,
        bwa_db = config["database_a"]["bwa_db"],
        reftype = config["database_a"]["type"]
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/map2ref.sh -t 16 -i {input.fwd} -n {input.rev} -r {params.bwa_db} -o {params.outpref} -c {params.reftype}
        """

# map metagenomes to catalog B
rule map2ref_b:
    input:
        fwd = lambda wildcards: samp2path[wildcards.sample][0],
        rev = lambda wildcards: samp2path[wildcards.sample][1]
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_B+"_total.tab",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_B+"_unique.tab"
    params:
        outpref = OUTPUT_DIR+"/mapping/{sample}/{sample}_" + REF_NAME_B,
        bwa_db = config["database_b"]["bwa_db"],
        reftype = config["database_b"]["type"]
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/map2ref.sh -t 16 -i {input.fwd} -n {input.rev} -r {params.bwa_db} -o {params.outpref} -c {params.reftype}
        """

# parse final output
rule parse_cov:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_A+"_total.tab", sample=samp2path.keys())
    output:
        OUTPUT_DIR+"/summary/bwa_cov-est.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-cov.py -i {params.outdir} -o {output}
        """

rule parse_expcov:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_A+"_total.tab", sample=samp2path.keys())
    output:
        OUTPUT_DIR+"/summary/bwa_cov-exp.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-expcov.py -i {params.outdir} -o {output}
        """

rule parse_counts:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_"+REF_NAME_A+"_{type}.tab", sample=samp2path.keys(), type=wildcards.type)
    output:
        OUTPUT_DIR+"/summary/bwa_counts-{type}.csv"
    params:
        outdir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-counts.py -i {params.outdir} -f {wildcards.type} -o {output}
        """
