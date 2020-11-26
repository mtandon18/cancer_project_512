import subprocess
import sys
import os
import HTSeq

def run(acc_list, index_base, gtf_file, paired):
    text = open(acc_list).read()
    lines = text.splitlines()
    output_dir = "SAM_Results"
    os.mkdir(output_dir)
    output_dir2 = "Count_Tables"
    os.mkdir(output_dir2)

    # Iterate through all samples to produce separate SAM files for all of them
    for line in lines:
        # Output file name, preserves SRA ID
        filename = line + ".sam"

        # Run hisat to produce alignments
        idx = "-x " + str(index_base)
        acc = "--sra-acc " + str(line)
        out = "-S " + filename
        align_summary = "--summary-file summary.txt" 
        metrics = "--met-file metrics.txt"
        # supress SAM records for reads that don't align, 
        cmd = ["hisat2 -p", idx, acc, out, "--no-unal", summary, metrics]
        subprocess.run(cmd)

        # 
        new_filename = line + "_counts.txt"
        redirect = f"> new_filename"
        nonunique = "--nonunique all"
        cmd = ["htseq-count", filename, gtf_file, redirect, "-n 4"]
        if (paired == "p"):
            subprocess.run(["samtools sort -n -o", filename, "-O sam -@ 4"])
            cmd += ["-r name"]
        subprocess.run(cmd)


        # Move files to new directories
        os.replace(filename, os.path.join(output_dir,filename))
        os.replace(new_filename, os.path.join(output_dir2,new_filename))

    

acc_list = sys.argv[1]
index_base = sys.argv[2]
gtf_file = sys.argv[3]
paired = sys.argv[4]
run(fastq_dir, index_base, gtf_file, paired)
