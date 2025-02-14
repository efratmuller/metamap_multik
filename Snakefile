# snakemake workflow for mapping metagenomes against 2 reference databases

# To run:
# cd /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/scripts/metamap_multik
# conda activate snakemake
# snakemake --use-conda -k -j 50 --cluster-config cluster.yml --cluster 'sbatch -A {cluster.project} -p {cluster.queue} --ntasks={cluster.nCPU} --mem={cluster.mem} -o {cluster.output} --time={cluster.time} -J {cluster.name}'

# Dry run:
# snakemake -n

import os

# read configuration yml
configfile: "config.yml"
INPUT_FILE = config["input"]
OUTPUT_DIR = config["output"]
REF_DBS = list(config["databases"].keys())
DB_PROPERTIES = config["databases"]

# verify exec permissions
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
        expand(OUTPUT_DIR+"/summary/bwa_{ref_db}_counts_total.csv", ref_db=REF_DBS),
        expand(OUTPUT_DIR+"/summary/bwa_{ref_db}_counts_unique_filtered.csv", ref_db=REF_DBS),
        expand(OUTPUT_DIR+"/summary/bwa_{ref_db}_cov_est.csv", ref_db=REF_DBS),
        expand(OUTPUT_DIR+"/summary/bwa_{ref_db}_cov_exp.csv", ref_db=REF_DBS),
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_multi_mapped_reads.txt", sample=samp2path.keys()),
        OUTPUT_DIR+"/summary/bwa_mapping_stats.tab"

# map metagenomes to catalog B (performed once per run acceesion)
rule map2ref:
    input:
        fwd = lambda wildcards: samp2path[wildcards.sample][0],
        rev = lambda wildcards: samp2path[wildcards.sample][1]
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_stats.tab",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_aligned_reads_uniq.txt.gz",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_unique_filtered.tab",
        OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_total.tab"
    params:
        outpref = OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}",
        bwa_db = lambda wildcards: DB_PROPERTIES[wildcards.ref_db]["bwa_db"],
        reftype = lambda wildcards: DB_PROPERTIES[wildcards.ref_db]["type"],
        cov_thresh = lambda wildcards: DB_PROPERTIES[wildcards.ref_db]["breadth_threshold"], 
        cov_exp_ratio_thresh = lambda wildcards: DB_PROPERTIES[wildcards.ref_db]["breadth_obs_exp_ratio"] 
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/map2ref.sh -t 16 -i {input.fwd} -n {input.rev} -r {params.bwa_db} -o {params.outpref} -c {params.reftype} -b {params.cov_thresh} -e {params.cov_exp_ratio_thresh}
        """

# calculate overlap of mapped reads
rule get_mapping_overlaps:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_aligned_reads_uniq.txt.gz", sample=wildcards.sample, ref_db=REF_DBS)
    output:
        OUTPUT_DIR+"/mapping/{sample}/{sample}_multi_mapped_reads.txt"
    shell:
        """
        zcat {input[0]} {input[1]} | sort | uniq -d > {output}
        """

# compute breadth of coverage and expected coverage per genome from each catalog
rule parse_cov:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_total.tab", sample=samp2path.keys(), ref_db=wildcards.ref_db)
    output:
        OUTPUT_DIR+"/summary/bwa_{ref_db}_cov_est.csv",
        OUTPUT_DIR+"/summary/bwa_{ref_db}_cov_exp.csv"
    params:
        in_dir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-cov.py -i {params.in_dir} -d {wildcards.ref_db} -o {output[0]}
        tools/parse_bwa-expcov.py -i {params.in_dir} -d {wildcards.ref_db} -o {output[1]}
        """

# organize final count tables (both unique and total, for both reference databases) 
rule parse_read_counts:
    input:
        lambda wildcards: expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_{type}.tab", sample=samp2path.keys(), ref_db=wildcards.ref_db, type=wildcards.type)
    output:
        OUTPUT_DIR+"/summary/bwa_{ref_db}_counts_{type}.csv"
    params:
        in_dir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-counts.py -i {params.in_dir} -f {wildcards.type} -d {wildcards.ref_db} -o {output}
        """

# collect statistics per sample into a single table
rule collect_stats:
    input:
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_multi_mapped_reads.txt", sample=samp2path.keys()),
        expand(OUTPUT_DIR+"/mapping/{sample}/{sample}_{ref_db}_stats.tab", sample=samp2path.keys(), ref_db=REF_DBS)
    output:
        OUTPUT_DIR+"/summary/bwa_mapping_stats.tab"
    params:
        in_dir = OUTPUT_DIR+"/mapping"
    conda:
        "envs/metamap.yml"
    shell:
        """
        tools/parse_bwa-stats.py -i {params.in_dir} -o {output}
        """