'''
This script takes a directory of fastq files as input and creates a corresponding
directory of SAM files once alignment is completed using HiSat. The naming convention
of the original files is retained. 
'''

import subprocess
import sys
import os

def run(folder):
    for filename in os.listdir(folder):
        result = filename.replace(".txt", "SAM.txt")
        subprocess.call(f"hisat2")




folder = sys.argv[1]
run(folder)
