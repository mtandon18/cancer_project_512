#!/usr/bin/env python3
import subprocess
import sys
import os
import HTSeq

def run(acc_list, index_base, gtf_file, paired, prefix):
    text = open(acc_list).read()
    lines = text.splitlines()
    output_dir = prefix + "_SAM_Results"
    os.mkdir(output_dir)
    output_dir2 = prefix + "_Count_Tables"
    os.mkdir(output_dir2)
    threads = "8"

    # Iterate through all samples to produce separate SAM files for all of them
    for line in lines:
        # Output file name, preserves SRA ID
        filename = line + ".sam"

        # Run hisat to produce alignments
        idx = "-x " + str(index_base)
        acc = "--sra-acc " + str(line)
        out = "-S " + filename
        summary = "--summary-file " + prefix + "_summary.txt"
        metrics = "--met-file " + prefix + "_metrics.txt"
        # supress SAM records for reads that don't align,
        cmd = ("hisat2 -p " + threads + " " + idx + " " + acc + " " + out + " --no-unal " + summary + " " + metrics)
        subprocess.run(cmd.split(" "))

        #
        new_filename = line + "_counts.txt"
        redirect = f"> {new_filename}"
        nonunique = "--nonunique all"
        cmd = "htseq-count " + filename + " " + gtf_file + " " + redirect + " -n " + threads
        if (paired == "p"):
            subprocess.run(("samtools sort -n -o " + filename + " -O sam -@ " + threads).split(" "))
            cmd += " -r name"
        subprocess.run(cmd.split(" "))


        # Move files to new directories
        os.replace(filename, os.path.join(output_dir,filename))
        os.replace(new_filename, os.path.join(output_dir2,new_filename))



acc_list = sys.argv[1]
index_base = sys.argv[2]
gtf_file = sys.argv[3]
paired = sys.argv[4]
prefix = sys.argv[5]
run(acc_list, index_base, gtf_file, paired, prefix)
