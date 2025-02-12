#!/usr/bin/env python

# Collect and merge all mapping statistics from all samples and all reference databases.

# To run:
# cd /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/scripts/metamap_multik
# tools/parse_bwa-stats.py -i /home/em2035/rds/hpc-work/mapping_tests/output_test/mapping -o /home/em2035/rds/hpc-work/mapping_tests/output_test/summary/all_stats_test.tsv

import os
import argparse
import pandas as pd

def parse_stats_filename(filename):
    """Extract run_accession and db_name from the filename."""
    parts = filename.rsplit("_", 2)  # Split from the right, expect 3 parts
    if len(parts) == 3 and parts[2] == "stats.tab":
        run_accession, db_name = parts[0], parts[1]
        return run_accession, db_name
    return None, None

def process_stats_files(input_directory):
    """Find, parse, and merge all _stats.tab files."""

    # Initialize a dictionary to store concatenated statistics per run_accession (across all reference DBs)
    stats_per_run = {}
    
    # Iterate over stats files
    for root, _, files in os.walk(input_directory):
        stat_dfs_per_run = {}

        # First, for each run accession, collect all statistics tables in a list
        for file in files:
            if file.endswith("_stats.tab"):
                print(f"Processing: {file}")
                run_accession, db_name = parse_stats_filename(file)
                if not run_accession or not db_name:
                    continue
                
                # Read stats
                df = pd.read_csv(os.path.join(root, file), sep="\t")

                # Add db_name column
                df["db_name"] = db_name  
                
                # Reorder 
                df = df[['category', 'db_name', 'key', 'value']]
                
                # Collect all files belonging to the same sample
                if run_accession not in stat_dfs_per_run:
                    stat_dfs_per_run[run_accession] = []
                stat_dfs_per_run[run_accession].append(df)
        
        # Next, concatenate all statistics tables for each run accession, to get one table with all statistics
        for run_accession, dfs in stat_dfs_per_run.items():
            concatenated_df = pd.concat(dfs, ignore_index=True)
            concatenated_df = concatenated_df.rename(columns={"value": run_accession}) # Rename value column
            stats_per_run[run_accession] = concatenated_df # Save
    
    # Merge all samples into a final table
    all_stats = None
    for stat_df in stats_per_run.values():
        if all_stats is None:
            all_stats = stat_df
        else:
            all_stats = all_stats.merge(stat_df, on=["category", "db_name", "key"], how="left")
    
    return(all_stats)

def count_overlapping_reads(input_directory):
    # Initialize a dictionary to store overlapping read counts
    overlap_reads_per_run = {}
    
    # Iterate over stats files
    for root, _, files in os.walk(input_directory):
        # First, for each run accession, collect all statistics tables in a list
        for filename in files:
            if filename.endswith("_multi_mapped_reads.txt"):
                print(f"Processing: {filename}")
                run_accession = filename.rsplit("_", 3)[0]  # Extract sample name

                # Count number of reads (lines)
                with open(os.path.join(root, filename), "r") as f:
                    line_count = sum(1 for _ in f)
                    overlap_reads_per_run[run_accession] = line_count

    # Convert to DataFrame
    df = pd.DataFrame([overlap_reads_per_run])
    df["key"] = "7_reads_mapped_to_both_catalogs"
    df["category"] = "READS"
    df["db_name"] = None
    return(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BWA results to extract counts')
    parser.add_argument('-i', dest='in_folder', help='Input folder (files should be in subdirectories)', required=True)
    parser.add_argument('-o', dest='out_file', help='Output TSV file', required=True)
    args = parser.parse_args()
    
    # Parse and collect main statistics
    all_stats = process_stats_files(args.in_folder)

    # Add statistics about reads mapped to both catalogs
    overlap_stats = count_overlapping_reads(args.in_folder)

    # Join the two tables
    overlap_stats = overlap_stats[all_stats.columns] # Reorder columns to match all_stats
    all_stats = pd.concat([all_stats, overlap_stats], ignore_index=True)

    if all_stats is not None:
        all_stats.to_csv(args.out_file, sep="\t", index=False)
        print(f"Merged file saved to {args.out_file}")
    else:
        print("No _stats.tab files found.")
