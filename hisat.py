#!/usr/bin/env python3
import subprocess
import sys
import os
import HTSeq

def run(fastq_files, index_base, gtf_file, paired, prefix):
    print(prefix, file=sys.stderr)
    output_dir = "SAM_Results"
    output_dir2 = "Count_Tables"
    threads = "8"

    # Output file name, preserves SRA ID
    filename = os.path.join(output_dir, prefix + ".sam")

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
    print("Running hisat2", file=sys.stderr)
    subprocess.run(cmd.split(" "))
    print("Ran hisat2, running htseq", file=sys.stderr)

    # Run htseq to produce count tables
    new_filename = os.path.join(output_dir2, prefix + "_counts.txt")
    new_file = open(new_filename, "w")
    nonunique = "--nonunique all"
    cmd = "htseq-count " + filename + " " + gtf_file + " -n " + threads
    filename_sorted = filename.replace(".sam", "_sort.sam")
    if (paired == "p"):
        subprocess.run(("samtools sort -n -o " + filename_sorted + " -O sam -@ " + threads + " " + filename).split(" "))
        os.remove(filename)
        filename = filename_sorted
        cmd += " -r name"
    subprocess.run(cmd.split(" "), stdout=new_file)
    print("Ran htseq", file=sys.stderr)
    new_file.close()

    # Move files to new directories
    #os.replace(filename_sorted, os.path.join(output_dir,filename_sorted))
    #os.replace(new_filename, os.path.join(output_dir2,new_filename))



fastq_files = sys.argv[1]
index_base = sys.argv[2]
gtf_file = sys.argv[3]
paired = sys.argv[4]
prefix = sys.argv[5]
run(fastq_files, index_base, gtf_file, paired, prefix)
