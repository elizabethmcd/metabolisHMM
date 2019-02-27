#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
from Bio import SeqIO, SearchIO

os.mkdir("out")
os.mkdir("results")
genomes=glob.glob("genomes/*.faa")
marker=sys.argv[1]
FNULL = open(os.devnull, 'w')
prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
dir=prot
os.mkdir("out/"+dir)

# Run HMM for a single marker
print("Searching for " + prot + " marker in genome set...")
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    outname= "out/"+dir+"/"+name + ".out"
    cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
    subprocess.call(cmd, stdout=FNULL)

# Parse HMM file 
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
print("Aligning hits...")
fastas = glob.glob("results/*.faa")
for fasta in fastas:
    outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
    output= "results/"+outname+".aln"
    musc_cmd = ["muscle","-quiet","-in",fasta,"-out",output]
    subprocess.call(musc_cmd)

# Make tree 
print("Constructing phylogeny...")
marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
alignment_file = "results/"+marker_name+".aln"
output_tree="results/"+marker_name+".tre"
tree_cmd = ["FastTree",alignment_file,">",output_tree]
subprocess.call(tree_cmd)