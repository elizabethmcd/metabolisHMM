#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
import pandas as pd 
from Bio import SearchIO

# Arguments 
DIR=sys.argv[1]

# Setup
genomes=glob.glob("genomes/*.faa")
markers=glob.glob(os.path.join(DIR, '*.hmm'))
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

# Parse HMM file to results matrix/dataframe
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

# Reformat dataframe in order of marker function, find markers in none of the genomes and input NaN for all
all_cols=[]
absent_cols=[]
existing_markers = df.columns
all_markers=glob.glob("metabolic_markers/*.hmm")
for marker in all_markers:
    prot=os.path.basename(marker).replace(".hmm", "")
    if prot not in existing_markers:
        all_cols.append(prot)
for col in existing_markers:
    all_cols.append(col)
df_all=df.reindex(columns=all_cols)
df_all.fillna(0, inplace=True)
df_all.to_csv("results/metabolic-marker-results.csv")