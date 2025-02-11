#!/usr/bin/env python

# Collect and merge all mapping statistics from all samples and all reference databases.

# To run:
# cd /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/scripts/metamap_multik
# tools/parse_bwa-stats.py -i /home/em2035/mapping_test/output_test/mapping -o /home/em2035/mapping_test/output_test/summary/all_stats_test.tsv

import os
import argparse
import pandas as pd

def parse_filename(filename):
    """Extract run_accession and db_name from the filename."""
    parts = filename.rsplit("_", 2)  # Split from the right, expect 3 parts
    if len(parts) == 3 and parts[2] == "stats.tab":
        run_accession, db_name = parts[0], parts[1]
        return run_accession, db_name
    return None, None

def process_stats_files(input_directory, output_file):
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
                run_accession, db_name = parse_filename(file)
                if not run_accession or not db_name:
                    continue
                
                # Read stats
                file_path = os.path.join(root, file)
                df = pd.read_csv(file_path, sep="\t")

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

            # Rename value column
            concatenated_df = concatenated_df.rename(columns={"value": run_accession})  

            # Save
            stats_per_run[run_accession] = concatenated_df
    
    # Merge all samples into a final table
    all_data = None
    for stat_df in stats_per_run.values():
        if all_data is None:
            all_data = stat_df
        else:
            all_data = all_data.merge(stat_df, on=["category", "db_name", "key"], how="left")
    
    if all_data is not None:
        all_data.to_csv(output_file, sep="\t", index=False)
        print(f"Merged file saved to {output_file}")
    else:
        print("No _stats.tab files found.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BWA results to extract counts')
    parser.add_argument('-i', dest='in_folder', help='Input folder (files should be in subdirectories)', required=True)
    parser.add_argument('-o', dest='out_file', help='Output TSV file', required=True)
    args = parser.parse_args()
    
    process_stats_files(args.in_folder, args.out_file)

