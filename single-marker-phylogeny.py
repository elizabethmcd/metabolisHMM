#! /usr/bin/env python 

import os, sys
import glob 
import argparse
import subprocess 
from Bio import SeqIO, SearchIO

parser = argparse.ArgumentParser(description = "Create phylogeny of single marker")
parser.add_argument('--genome_dir', metavar='GENOMEDIR', help='Directory where genomes to be screened are held')
parser.add_argument('--marker', metavar='MARKER', help="Location of single marker to run analysis on")
parser.add_argument('--phylogeny', metavar='PHY', help="fastree or raxml, choose one")
parser.add_argument("--threads",metavar='THREADS',help="number of threads for tree making")

args = parser.parse_args()
GENOMEDIR = args.genome_dir
MARKER = args.marker
PHYTOOL = args.phylogeny
THREADS = args.threads

os.mkdir("out")
os.mkdir("results")
genomes=glob.glob(os.path.join(GENOMEDIR, '*.faa'))
marker=MARKER
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
if PHYTOOL == 'fastree':
    print("Calculating tree using FastTree...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    alignment_file = "results/"+marker_name+".aln"
    output_tree="results/"+marker_name+".tre"
    tree_cmd = ["FastTree","-out",output_tree,alignment_file]
    subprocess.call(tree_cmd)
elif PHYTOOL == "raxml":
    print("Calculating tree with RaxML... be patient...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    outname= marker_name+"raxml"
    fileIn="results/"+marker_name+".aln"
    raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
    os.system(raxCmd)