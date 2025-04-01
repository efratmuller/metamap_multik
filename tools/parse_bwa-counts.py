#!/usr/bin/env python

# To test:
# cd /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/scripts/metamap_multik
# tools/parse_bwa-counts.py -f "unique_filtered" -d "uhgv" -i /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/mapping/positive_controls_test/mapping -o /rds/project/rds-aFEMMKDjWlo/emuller/phages_project/mapping/positive_controls_test/summary/bwa_uhgv_counts_unique_filtered.csv

import glob
import sys
import os
import argparse
import pprint
import csv

def collectCounts(folder, ftype, db_name, out_file):
    counts = {}
    file_count = 0
    print("Collecting read counts from the following files:")
    for f in glob.glob(os.path.join(folder, f"*/*_{db_name}_{ftype}.tab")):
        file_count += 1
        print("%i\t%s" % (file_count, f))
        with open(f, "r") as bwa_data:
            linen = 0
            for line in bwa_data:
                linen += 1
                cols = line.strip("\n").split("\t")

                # Record sample ID (extract it from the header line)
                if linen == 1:
                    name = os.path.basename(cols[2])
                    name = name.replace("_" + db_name + "_Counts", "")
                    counts[name] = {}
                # For non-header rows, simply record read counts
                else:
                    counts[name][cols[0]] = cols[2]

    # Now convert the dictionary into a table:
    ## Collect all unique genomes
    genomes = set()
    for sample in counts.values():
        genomes.update(sample.keys())

    ## Prepare output
    header = ['Genome'] + list(counts.keys())
    rows = []

    ## Step 3: Populate rows with genome names and values
    for genome in genomes:
        row = [genome]
        for sample in counts:
            row.append(counts[sample].get(genome, 0))  # Use 0 if the genome is not in the sample
        rows.append(row)

    # Step 4: Write to a CSV file
    with open(out_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)  # Write the header
        writer.writerows(rows)  # Write the rows

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BWA results to extract counts')
    parser.add_argument('-i', dest='in_folder', help='Input folder (files should be in subdirectories)', required=True)
    parser.add_argument('-f', dest='ftype', help='\'unique\' or \'unique_filtered\' or \'total\' counts', required=True)
    parser.add_argument('-d', dest='db_name', help='Name of DB of reference genomes', required=True)
    parser.add_argument('-o', dest='out_file', help='Output CSV file', required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    print("User-provided arguments:")
    pprint.pprint(vars(args))
    collectCounts(args.in_folder, args.ftype, args.db_name, args.out_file)

