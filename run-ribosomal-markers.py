#! /usr/bin/env python 

import glob
import subprocess 
import os
import pandas as pd 
from Bio import SeqIO, SearchIO
from Bio.Nexus import Nexus

# Setup
genomes = glob.glob("genomes/*.faa")
ribos = glob.glob("ribosomal_markers/*.hmm")
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

# setup hmmsearch run
for genome in genomes:
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir("out/"+dir)
    for ribo in ribos: 
        prot=os.path.basename(ribo).replace(".hmm", "").replace("_bact", "").strip().splitlines()[0]
        outname= "out/"+dir+"/"+name + "-" + prot + ".out"
        cmd = ["hmmsearch", "--tblout="+outname, ribo, genome]
        subprocess.call(cmd, stdout=FNULL)
        print("Running HMMsearch on " + name + " and " + prot + " marker")

# parsing the hmm output file and get into results file

result_dirs = os.walk("out/")
output = ("results/all-results-loci.txt")
with open(output, "a") as f:
    for path, dirs, files in result_dirs:
        for file in files:
            genome = file.split("-")[0]
            prot = file.replace(".out", "").split("-")[1]
            result = "out/"+genome+"/"+file
            with open(result, "rU") as input:
                for qresult in SearchIO.parse(input, "hmmer3-tab"):
                    hits = qresult.hits
                    num_hits = len(hits)
                    if num_hits > 0:
                        for i in range(0, 1): # takes only the first, best hit
                            hit_id = hits[i].id
                            f.write(genome + "\t" + prot + "\t" + hit_id + "\n")

# Parse locus tags from results file, get sequences 
loci_list = pd.read_table("results/all-results-loci.txt", sep="\t", names=['genome', 'marker', 'locus'])



# Make alignment file
