#! /usr/bin/env python 

import glob
import subprocess 
import os, sys
import pandas as pd 
from Bio import SeqIO, SearchIO

# Setup
# genomes = glob.glob("genomes/*.faa")
# ribos = glob.glob("ribosomal_markers/*.hmm")
# os.mkdir("out")
# os.mkdir("results")
# FNULL = open(os.devnull, 'w')

# # setup hmmsearch run
# for genome in genomes:
#     name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
#     dir=name
#     os.mkdir("out/"+dir)
#     for ribo in ribos: 
#         prot=os.path.basename(ribo).replace(".hmm", "").replace("_bact", "").strip().splitlines()[0]
#         outname= "out/"+dir+"/"+name + "-" + prot + ".out"
#         cmd = ["hmmsearch", "--tblout="+outname, ribo, genome]
#         subprocess.call(cmd, stdout=FNULL)
#         print("Running HMMsearch on " + name + " and " + prot + " marker")

# Parse HMM outputs
prot_list=['rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8', 'rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8']
result_dirs = os.walk("out/")
for prot in prot_list:
    for path, dirs, files in result_dirs: 
        for file in files:
            genome=file.split("-")[0]
            marker=file.replace(".out", "").split("-")[1]
            result="out/"+genome+"/"+file
            output="results/"+marker+".faa"
            genome_file = "genomes/"+genome+".faa"
            with open(output, "a") as outf:
                with open(genome_file, "r") as input_fasta:
                    with open(result, "rU") as input:
                        for qresult in SearchIO.parse(input, "hmmer3-tab"):
                            hits=qresult.hits
                            num_hits=len(hits)
                            if num_hits >0:
                                for i in range(0, 1):
                                    hit_id=hits[i].id
                            for record in SeqIO.parse(input_fasta, "fasta"):
                                if record.id in hit_id:
                                    outf.write(">"+genome+"_"+hit_id+"\n"+str(record.seq)+"\n")

# Make alignments for each marker
