#!/usr/bin/env python3
import subprocess
import sys
import os
import HTSeq

def run(fastq_files, index_base, gtf_file, paired, prefix):
    output_dir = "SAM_Results"
    output_dir2 = "Count_Tables"
    threads = "8"

    # Output file name, preserves SRA ID
    filename = prefix + ".sam"

    # Run hisat to produce alignments
    fastq = None
    if (paired == "p"):
        one, two = fastq_files.split(",")
        fastq = f"-1 {one} -2 {two}"
    else:
        fastq = f"-U {fastq_files}"

    idx = "-x " + str(index_base)
    out = "-S " + filename
    summary = "--summary-file " + prefix + "_summary.txt"
    metrics = "--met-file " + prefix + "_metrics.txt"
    # supress SAM records for reads that don't align,
    cmd = ("hisat2 -p " + threads + " " + idx + " " + fastq + " " + out + " --no-unal " + summary + " " + metrics)
    subprocess.run(cmd.split(" "))

    # Run htseq to produce count tables
    new_filename = prefix + "_counts.txt"
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



fastq_files = sys.argv[1]
index_base = sys.argv[2]
gtf_file = sys.argv[3]
paired = sys.argv[4]
prefix = sys.argv[5]
run(fastq_files, index_base, gtf_file, paired, prefix)
