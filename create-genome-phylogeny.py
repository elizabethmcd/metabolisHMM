#! /usr/bin/env python 

import glob
import argparse
import subprocess 
import os, sys
import pandas as pd 
from Bio import SeqIO, SearchIO

# Arguments

parser = argparse.ArgumentParser(description = "Create full archaeal/bacterial genome phylogenies based off specific ribosomal protein markers")
parser.add_argument('--genome_dir', metavar='GENOMEDIR', help='Directory where genomes to be screened are held')
parser.add_argument('--domain', metavar='DOMAIN', help='archaea, bacteria (choose one, need to be assessed separately)')

args = parser.parse_args()
GENOMEDIR = args.genome_dir
GENOMEFILES = GENOMEDIR + "/*.faa"
DOMAIN = args.domain

# Setup
genomes = glob.glob(GENOMEDIR)
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

bacteria_list = ['rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8', 'rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8']
archaea_list = ['rpL14', 'rpL15', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS17', 'rpS19', 'rpS3', 'rpS8', 'rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8']

# setup hmmsearch run for archaea or bacteria
for genome in GENOMFILES:
    name=os.path.basename().replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir("out/"+dir)
    if DOMAIN == 'archaea':
        for prot in archaea_list:
            marker ="ribosomal_markers/"+prot+"_bact.hmm"
            outname= "out/"+dir+"/"+name + "-" + prot + ".out"
            cmd = ["hmmsearch", "--tblout="+outname, marker, genome]
            subprocess.call(cmd, stdout=FNULL)
            print("Running HMMsearch on " + name + " and " + prot + " marker")
    elif DOMAIN == 'bacteria':
        for prot in bacteria_list:
            marker="ribosomal_markers/"+prot+"_bact.hmm"
            outname= "out/"+dir+"/"+name + "-" + prot + ".out"
            cmd = ["hmmsearch", "--tblout="+outname, marker, genome]
            subprocess.call(cmd, stdout=FNULL)
            print("Running HMMsearch on " + name + " and " + prot + " marker")
    
# Parse HMM outputs
print("Parsing results...")
if DOMAIN == 'archaea':    
    prot_list=['rpL14', 'rpL15', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6','rpS17', 'rpS19', 'rpS3', 'rpS8', 'rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8']
else:
    prot_list=['rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8', 'rpL14', 'rpL15', 'rpL16', 'rpL18', 'rpL2', 'rpL22', 'rpL24', 'rpL3', 'rpL4', 'rpL5', 'rpL6', 'rpS10', 'rpS17', 'rpS19', 'rpS3', 'rpS8']
result_dirs = os.walk("out/")
for prot in prot_list:
    for path, dirs, files in result_dirs: 
        for file in files:
            genome=file.split("-")[0]
            marker=file.replace(".out", "").split("-")[1]
            result="out/"+genome+"/"+file
            output="results/"+marker+".faa"
            genome_file = GENOMEDIR+genome+".faa"
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
                                    outf.write(">"+genome+"\n"+str(record.seq)+"\n")
                       
# Make alignment file
print("Aligning concatenated hits...")
fastas = glob.glob("results/*.faa")
for fasta in fastas:
    outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
    output= "results/"+outname+".aln"
    musc_cmd = ["muscle","-quiet","-in",fasta,"-out",output]
    subprocess.call(musc_cmd)

# Concatenate alignments 
