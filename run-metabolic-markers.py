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
print("Parsing all results...")
all_dicts={}
result_dirs = os.walk("out/")
for path, dirs, files in result_dirs:
    genome_dict={}
    for file in files:
        genome = file.split("-")[0]
        prot = file.replace(".out", "").split("-")[1]
        result = "out/"+genome+"/"+file
        with open(result, "rU") as input:
            for qresult in SearchIO.parse(input, "hmmer3-tab"):
                hits = qresult.hits
                num_hits = len(hits)
                genome_dict[prot] = num_hits
                all_dicts[os.path.basename(file).split("-")[0]]=genome_dict
df=pd.DataFrame.from_dict(all_dicts, orient="index", dtype=None)
df.to_csv("results/metabolic-markers-results.csv")
df.to_csv("results/metabolic-results.tsv", sep="\t")
