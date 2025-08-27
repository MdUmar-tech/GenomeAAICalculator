## AAI Calculator (Average Aminoacid Identity)
A Python tool for calculating Average Amino Acid Identity (AAI) between two or more genomes using reciprocal best BLAST hits.

## Description
This tool calculates the Average Amino Acid Identity (AAI) between genomic datasets by performing pairwise comparisons of protein sequences. It uses BLAST+ to find reciprocal best hits between genomes and computes AAI values based on the percentage identity of these matches.

## Features
Automated BLAST database creation for protein sequences

Reciprocal BLASTp searches between genome pairs

Calculation of AAI values with standard deviation

Multi-threading support for faster processing

CSV output of results

Option to clean up intermediate files

## Requirements
Python 3.6+

BLAST+ installed in your system PATH

pandas library

Installation
Ensure BLAST+ is installed on your system:

bash
# On Ubuntu/Debian
sudo apt-get install ncbi-blast+

# On macOS with Homebrew
brew install blast
Install required Python packages:

bash
pip install pandas
Clone or download this repository.

Usage
bash
python AAI_calculate.py -i input_directory [-o output_file] [-t threads] [--clean]
Parameters
-i, --input: Input directory containing FASTA files with protein sequences (required)

-o, --output: Output file name (default: 'output_final.csv')

-t, --threads: Number of CPU threads to use for BLAST (default: 1)

--clean: Remove intermediate BLAST database and output files

Example
bash
python AAI_calculate.py -i genomes/ -o aai_results.csv -t 4 --clean
Input Format
Place all protein sequence files in FASTA format (.faa extension) in the input directory. Each file should represent one genome.

Example file structure:

text
genomes/
    genome1.faa
    genome2.faa
    genome3.faa
Output
The tool generates several output files:

Pairwise RBBH files: CSV files containing reciprocal best hit data for each genome pair

Pairwise AAI results: CSV files with AAI values for each genome pair

Final output matrix: A comprehensive CSV file (default: output_final.csv) with all pairwise AAI values

The final output includes:

Genome identifiers

AAI percentage

Number of proteins used in calculation

Standard deviation of identity values

Algorithm
Create BLAST databases for each protein file

Perform pairwise reciprocal BLASTp searches between all genomes

Identify reciprocal best hits (RBBH)

Calculate AAI as the average identity of RBBH

Compile results into a matrix format

Author
Md Umar

Contact: arc.umar@cusat.ac.in

References
This implementation is based on the concept of Average Amino Acid Identity calculation as described in microbial genomics literature.

Notes
Ensure all input files have the .faa extension

The process may take significant time for large genomes or many comparisons

Using more threads (-t option) can speed up BLAST searches

The --clean option helps save disk space by removing intermediate files

License
This project is provided as-is without specific license information. Please contact the author for usage permissions.
