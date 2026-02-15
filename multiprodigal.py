import os
import subprocess

# Define input and output directories
input_dir = "all_fasta"
output_dir = "output"

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# List all files in the input directory
input_files = os.listdir(input_dir)

# Iterate through each input file
for file in input_files:
    if file.endswith(".fna"):
        # Formulate the command
        command = f"prodigal -i {os.path.join(input_dir, file)} -o {os.path.join(output_dir, file.split('.')[0] + '.genes')} -a {os.path.join(output_dir, file.split('.')[0] + '.proteins.faa')}"
        
        # Execute the command
        subprocess.run(command, shell=True)
