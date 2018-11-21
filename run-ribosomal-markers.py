#! /usr/bin/env python 

import glob
import subprocess 
import os
from Bio import SeqIO, SearchIO
from Bio.Nexus import Nexus

# Setup
genomes = glob.glob("genomes/*.faa")
ribos = glob.glob("ribosomal_markers/*.hmm")
os.mkdir("out")

# setup hmmsearch run
for genome in genomes:
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir("out/"+dir)
    for ribo in ribos: 
        prot=os.path.basename(ribo).replace(".hmm", "").replace("_bact", "").strip().splitlines()[0]
        outname= "out/"+dir+"/"+name + "-" + prot + ".out"
        print(outname)
        cmd = ["hmmsearch", "--domtblout="+outname, ribo, genome]
        subprocess.call(cmd)
        print("Running HMMsearch on " + name + " and " + prot + " marker")

# parsing the hmm output file and get into results file

result_dirs = os.walk("out/")
for path, dirs, files in result_dirs:
    for dir in dirs:
        genome=dir
        output=(genome+"-loci.txt")
    for file in files:
        prot=file.replace(".out", "")
        result= "out/" + genome + "/" + file
        with open(result, "rU") as input:
            for qresult in SearchIO.parse(input, "hmmsearch3-domtab"):
                hits = qresult.hits
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0, num_hits):
                        hit_id = hits[i].id
                        print(hit_id)


# Parse locus tags from results file, get sequences 

# Make alignment file
