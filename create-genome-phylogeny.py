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
parser.add_argument('--phylogeny', metavar='PHY', help="fastree or raxml, choose one")
parser.add_argument("--threads",metavar='THREADS',help="number of threads for tree making")

args = parser.parse_args()
GENOMEDIR = args.genome_dir
GENOMEFILES = GENOMEDIR + "/*.faa"
DOMAIN = args.domain
PHYTOOL = args.phylogeny
THREADS = args.threads

# Setup
genomes = glob.glob(GENOMEFILES)
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

bacteria_list = ['rpL14','rpL15','rpL16','rpL18','rpL22','rpL24','rpL2','rpL3','rpL4','rpL5','rpL6','rpS10','rpS17','rpS19','rpS3','rpS8']
archaea_list = ['rpL14','rpL15','rpL18','rpL22','rpL24','rpL2','rpL3','rpL4','rpL5','rpL6','rpS17','rpS19','rpS3','rpS8']

# setup hmmsearch run for archaea or bacteria
if DOMAIN == 'archaea':
    for genome in genomes:
        name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
        dir=name
        os.mkdir("out/"+dir)
        for prot in archaea_list:
            marker ="ribosomal_markers/"+prot+"_bact.hmm"
            outname= "out/"+dir+"/"+name + "-" + prot + ".out"
            cmd = ["hmmsearch", "--tblout="+outname, marker, genome]
            subprocess.call(cmd, stdout=FNULL)
            print("Running HMMsearch on " + name + " and " + prot + " marker")
elif DOMAIN == 'bacteria':
    for genome in genomes:
        name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
        dir=name
        os.mkdir("out/"+dir)
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
                    with open(result, "r") as input:
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
    mafft_cmd = "mafft --quiet "
    mafft_cmd += fasta+" > "+output
    os.system(mafft_cmd)

# Concatenate alignments
print("Concatenating alignments...")
alignments="results/*.aln"
outname="results/"+DOMAIN+"-ribo-concatenated-phylogeny.fasta"
cat_cmd = "catfasta2phyml.pl -f --concatenate "
cat_cmd += alignments+" > "+outname
os.system(cat_cmd)

# Create tree
if PHYTOOL == 'fastree':
    print("Calculating tree using FastTree...")
    fileIn="results/"+DOMAIN+"-ribo-concatenated-phylogeny.fasta"
    outname = "results/"+DOMAIN+"-fastTree-ribosomal-tree.tre"
    fastCmd = "FastTree "+fileIn+" > "+outname
    os.system(fastCmd)
elif PHYTOOL == "raxml":
    print("Calculating tree with RaxML... be patient...")
    outname= DOMAIN+"-raxml-ribo"
    outDir = "results/"
    fileIn="results/"+DOMAIN+"-ribo-concatenated-phylogeny.fasta"
    raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
    os.system(raxCmd)