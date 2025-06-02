#!/usr/bin/env python

# This script takes the matrix of all observed genome coverages (rows = genomes, columns = samples), 
#  an equivalant matrix of the expected coverage / observed coverage ratio, and the threshold used 
#  for each metric to infer presence/absence, and generate a matrix with the same dimensions that 
#  simply holds 1 for a genome inferred as present in a sample and 0 otherwise.

# To run:
# cd /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/scripts/metamap_multik
# tools/get_genome_presence_absence_matrices.py --csv1 ../../mapping/all_01042025/summary/bwa_uhgg_cov_est.csv --csv2 ../../mapping/all_01042025/summary/bwa_uhgg_cov_exp_obs_ratio.csv --t1 5 --t2 0.3 -o ../../mapping/all_01042025/summary/bwa_uhgg_presence_absence_test.csv

import pandas as pd
import numpy as np
import argparse

def process_csvs(csv1_path, csv2_path, threshold1, threshold2, output_path, chunksize=200):
    reader1 = pd.read_csv(csv1_path, chunksize=chunksize)
    reader2 = pd.read_csv(csv2_path, chunksize=chunksize)

    with open(output_path, 'w') as out_f:
        header_written = False
        i = 0
        # Read csvs, chunk by chunk
        for chunk1, chunk2 in zip(reader1, reader2):
            i += 1
            print(f"Processing chunk {i}")

            # Write the header line once
            if not header_written:

                # Verify that indeed the headers of the two csvs match and if not throw an error
                if len(chunk1.columns) != len(chunk2.columns):
                    raise ValueError("CSV files have a different number of columns.")
                if not all(chunk1.columns == chunk2.columns):
                    raise ValueError("CSV headers do not match.")
                out_f.write(','.join(chunk1.columns) + '\n')
                header_written = True

            # Also verify that the values in the 1st column (genome IDs) are the same
            if not (chunk1.iloc[:, 0].equals(chunk2.iloc[:, 0])):
                raise ValueError("First columns (keys) do not match in a chunk.")

            ids = chunk1.iloc[:, 0] # Extract genome IDs
            data1 = chunk1.iloc[:, 1:].astype(float)
            data2 = chunk2.iloc[:, 1:].astype(float)

            result = ((data1 >= threshold1) & (data2 >= threshold2)).astype(int)
            result.insert(0, chunk1.columns[0], ids) # Place the genome IDs in the leftmost column

            result.to_csv(out_f, header=False, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a presence (1) absence (0) matrix based on two large matrices with the same rows (genomes) and columns (samples), and two thresholds. Presence in the output matrix is set to 1 only if the corresponding cell in the 1st matrix is above threshold1, and the corresponding cell in the 2nd matrix is above threshold2.")
    parser.add_argument("--csv1", required=True, help="Path to the first CSV file")
    parser.add_argument("--csv2", required=True, help="Path to the second CSV file")
    parser.add_argument("--t1", type=float, required=True, help="Threshold for CSV1 values")
    parser.add_argument("--t2", type=float, required=True, help="Threshold for CSV2 values")
    parser.add_argument("--output", "-o", required=True, help="Path to the output CSV file")

    args = parser.parse_args()

    process_csvs(args.csv1, args.csv2, args.t1, args.t2, args.output)
    print("Done")
