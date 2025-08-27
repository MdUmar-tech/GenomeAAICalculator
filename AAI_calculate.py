#!/usr/bin/python3
"""
Calculate the Average Aminoacid Identity (POCP) between two or
more genomes 
The program was written based on ()
# Required:
BLAST+ installed in $PATH
# Usage:
$ python AAI_calculate.py -i input_dir [-o output_matrix.tab] [-t 8] [--clean]
# Options:
-i: input directory contained two or more than 2 translated genome files (suffix: .faa)
-o: output AAI_final.csv
-n: number of threads (optional, default: 3)
--clean: blast output and databases created by this program will be removed (optional)

"""
import csv
import glob
import itertools
import os
import subprocess
import pandas as pd
import argparse

__author__ = "Md Umar"
__contact__ = "arc.umar@cusat.ac.in"

parser = argparse.ArgumentParser(description='Calculates the Average Amino Acid Identity within a pair of genomes')
parser.add_argument('-i', '--input', dest='input_directory', type=str, required=True, help='path to the input directory that contains the RBBH files')
parser.add_argument('-o', '--output', dest='output_directory', type=str, default='output_final.csv', help='path to the output directory for the AAI results')

parser.add_argument('-t', '--threads', metavar='num_cpu', type=int, default=1, help='The number of CPU threads to employ for BLAST+, Default is 1')
parser.add_argument('--clean', metavar='clean_blast_db_output', dest='c', nargs="?", const=True, default=False, help='Delete unnecessary files made by this programme if added')

args = parser.parse_args()


#print("Running one-way comparisons." if not o["q"] else "")

import glob
import itertools
import os
import subprocess
import pandas as pd


def run_mkblastdb(fi, fo):
    cmd_para = [
                'makeblastdb',
                '-in', fi,
                '-dbtype', 'prot',
                '-parse_seqids',
                '-out', fo
                ]
    try:
        run = subprocess.call(cmd_para, stdout=subprocess.PIPE)
    except Exception as e:
        raise e


def run_blastp(q, db, o, n):
    cmd_para = [
                'blastp',
                '-query', q,
                '-out', o,
                '-db', db,
                '-num_threads', str(n),
                '-outfmt', "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                '-max_target_seqs', '1',
                 '-max_hsps', '1'
                ]
    print("Executing BLAST command:")
    print(" ".join(cmd_para))
    try:
        run = subprocess.call(cmd_para, stdout=subprocess.PIPE)
    except Exception as e:
        raise e


def get_rbbh(genome_pair):
    blastout_name1 = genome_pair[0] + '--' + os.path.basename(genome_pair[1]) + '.AAIout'
    blastout_name2 = genome_pair[1] + '--' + os.path.basename(genome_pair[0]) + '.AAIout'
    run_blastp(genome_pair[0], genome_pair[1]+'_AAI', blastout_name1, args.threads)
    run_blastp(genome_pair[1], genome_pair[0]+'_AAI', blastout_name2, args.threads)
    fwd_results = pd.read_csv(blastout_name1, sep="\t", header=None)
    rev_results = pd.read_csv(blastout_name2, sep="\t", header=None)
    print(f"First 5 lines of {blastout_name1}:")
    print(fwd_results.head())

    print(f"First 5 lines of {blastout_name2}:")
    print(rev_results.head())

    headers = ["query" ,"subject", "identity", "length", "mismatch", "gapopen", "qstart",
               "qend", "sstart", "send", "evalue", "bitscore"]
    fwd_results.columns = headers
    rev_results.columns = headers
    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
            left_on='subject', right_on='query',
            how='outer')

    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
    return rbbh


def write_csv(genome_pair, rbbh):
    genome1 = os.path.basename(genome_pair[0]).split(".")[0]
    genome2 = os.path.basename(genome_pair[1]).split(".")[0]
    outfilename = f"{genome1}_{genome2}_RBBH.csv"
    rbbh.to_csv(outfilename, index=True)


def calculate_aai(rbbh_file):
    if not os.path.isfile(rbbh_file):
        print(f"RBBH file {rbbh_file} not found. Skipping AAI calculation for rbbh.")
        return None

    parts = rbbh_file.split(".faa_")

    if len(parts) != 3:
        raise ValueError("Input file name does not follow the expected format")

    genome1 = parts[0]
    genome2 = parts[1].split(".")[0]

    if not os.path.isfile(rbbh_file):
        print(f"RBBH file {rbbh_file} not found. Skipping AAI calculation for rbbh.")
        return None

    reciprocal_hits = pd.read_csv(rbbh_file, sep=",")
    n_protein = len(reciprocal_hits)
    if n_protein == 0:
        print(f"No reciprocal best hits found in {rbbh_file}. Skipping AAI calculation.")
        return None

    aai = float(reciprocal_hits["identity"].sum() / n_protein)
    aai = round(aai, 2)  # Format to two decimal places
    sd = float(reciprocal_hits["identity"].std())
    sd = round(sd, 2)  # Format to two decimal places

    # Define the header for the CSV file
    header = ["genome1", "genome2", "AAI", "n_protein", "sd"]

    # Define the data to be written to the CSV file
    data = [[genome1, genome2, aai, n_protein, sd]]

    outfilename = f"{genome1}_{genome2}_AAI_results.csv"

    # Open the CSV file in append mode and write the header and data to a new row
    with open(outfilename, "a", newline="") as f:
        writer = csv.writer(f)
        if os.path.getsize(outfilename) == 0:
            writer.writerow(header)  # Write the header to the CSV file
        writer.writerows(data)  # Write the data to the CSV file

        
def write_final_csv(data, header):
    # Open the CSV file in append mode and write the header and data to a new row
    with open("output_final.csv", "a", newline="") as f:
        writer = csv.writer(f)
        # Check if the file is empty
        if f.tell() == 0:
            writer.writerow(header)  # Write the header to the CSV file
        writer.writerows(data)  # Write the data to the CSV file

    

def main():
    print("Starting the main function...")
    genomes = glob.glob(os.path.join(args.input_directory,'*.faa'))
    print(f"Found genomes: {genomes}")
    for genome in genomes:
        run_mkblastdb(genome, genome+'_AAI')
        print(f"BLAST database created for {genome}")

    genome_pairs = list(itertools.combinations(genomes, 2))
    print(f"Found genomes: {genomes}")
    for genome_pair in genome_pairs:
        print(f"Processing genome pair: {genome_pair}")
        rbbh = get_rbbh(genome_pair)
        print(f"RBBH results for {genome_pair}:")
        print(rbbh.head())
        write_csv(genome_pair, rbbh)
        print(f"CSV written for {genome_pair}")
    # Initialize a list to store the processed file names
    processed_files = []

    # Load the list of processed files from a file, if it exists
    if os.path.isfile("processed_files.txt"):
        with open("processed_files.txt", "r") as f:
            processed_files = f.read().splitlines()

    # Process multiple CSV files using the glob module
    inputfiles = glob.glob("*.csv")

    for inputfile in inputfiles:
        if ".faa_" in inputfile and inputfile not in processed_files:
            calculate_aai(inputfile)
            # Add the processed file name to the list
            processed_files.append(inputfile)
            with open("processed_files.txt", "w") as f:
                f.write("\n".join(processed_files))
            # Write the list of processed files to
    
    #inputfiles = [f"{genome_pair[0]}_{genome_pair[1]}__RBBH.csv" for genome_pair in genome_pairs]
    
    #write_final_csv(inputfiles, outputfile)
    #write_final_csv(outputfile)


if __name__ == '__main__':
    main()
