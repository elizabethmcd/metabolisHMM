#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
import pandas as pd 
from Bio import SearchIO

genomes=glob.glob("genomes/*.faa")
marker=sys.argv[1]
FNULL = open(os.devnull, 'w')

# Run HMM for a single marker
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    dir=prot
    os.mkdir("out/"+dir)
    outname= "out/"+dir+"/"+name + ".out"
    cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
    subprocess.call(cmd, stdout=FNULL)
    print("Searching for " + prot + " marker in genome set")

# Parse HMM file to results matrix/dataframe
print("Parsing all results...")
result_dir = os.walk("out/"+dir)
for path, dirs, files in result_dir:
    for file in files:
        genome = file.replace(".out", "").strip().splitlines()[0]
        result = "out/"+dir+"/"+file
        output="results/"+dir+".faa"
        genome_file="genomes/"+genome+".faa"
        with open(output, "a") as outf:
            with open(genome_file, "r") as input_fasta:
                with open(result, "rU") as input:
                    for qresult in SearchIO.parse(input, "hmmer3-tab"):
                        hits = qresult.hits
                        num_hits = len(hits)
                        if num_hits>0:
                            for i in range(0,1):
                                hit_id=hits[i].id
                            for record in SeqIO.parse(input_fasta, "fasta"):
                                if record.id in hit_id:
                                    outf.write(">"+genome+"\n"+str(record.seq)+"\n")

# Align hits 


# Make tree 