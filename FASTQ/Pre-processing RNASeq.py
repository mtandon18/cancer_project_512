import subprocess
import sys
import os

def run(acc_list, index_base):
    text = open(acc_list).read()
    lines = text.splitlines()
    output_dir = "SAM_Results"
    os.mkdir(output_dir)

    # Iterate through all samples to produce separate SAM files for all of them
    for line in lines:
        # Output file name, preserves SRA ID
        filename = line + ".sam"

        # Run command in shell
        cmd = ["hisat2 -x ", index_base, "--sra-acc ", line, "-S ", filename]
        subprocess.run(cmd)
        
        # Move file to new directory
        os.replace(filename, os.path.join(output_dir,filename))

acc_list = sys.argv[1]
index_base = sys.argv[2]
run(fastq_dir, index_bas)
