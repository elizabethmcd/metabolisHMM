#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
import pandas as pd 
from Bio import SearchIO

# Setup
genomes=glob.glob("genomes/*.faa")
markers=glob.glob("metabolic_markers/*.hmm")
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

# Run HMMs
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir("out/"+dir)
    for marker in markers:
        prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        outname= "out/"+dir+"/"+name + "-" + prot + ".out"
        cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
        subprocess.call(cmd, stdout=FNULL)
        print("Running HMMsearch on " + name + " and " + prot + " marker")

# Parse HMM file to results matrix 
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
                        f.write(genome + "\t" + prot + "\t" + str(num_hits) + "\n")
