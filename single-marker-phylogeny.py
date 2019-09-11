#! /usr/bin/env python 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# single-marker-phylogeny : to create a phylogeny of a single HMM marker
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import os, sys, glob, subprocess, argparse
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO, SeqIO

# Arguments and Directory setup
parser = argparse.ArgumentParser(description = "Create phylogeny of single marker")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
required.add_argument('--input', metavar='INPUT', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help='Directory to store results and intermediate files')
required.add_argument('--marker', metavar='MARKER', help="Location of single marker to run analysis on")
optional.add_argument('--list', default='hit_list.txt', help="Output list of hits locus tags")
required.add_argument('--phylogeny', metavar='PHY', help="fastree or raxml, choose one")
optional.add_argument("--threads",metavar='THREADS',help="number of threads for tree making")

# if no arguments given, print help message
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# version to print
def version():
    versionFile = open('VERSION')
    return versionFile.read()
VERSION = version()

args = parser.parse_args()
GENOMEDIR = args.input
OUTPUT = args.output
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
MARKER = args.marker
PHYTOOL = args.phylogeny
THREADS = args.threads

os.makedirs(out_intm)
os.makedirs(out_results)
genomes=glob.glob(os.path.join(GENOMEDIR, '*.faa'))
marker=MARKER
FNULL = open(os.devnull, 'w')
prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
dir=prot
os.makedirs(OUTPUT + "/out/"+dir)

# optional write out hit list
OUT_LIST_PATH = OUTPUT + "/results/" + args.list
OUT_LIST = open(OUT_LIST_PATH, "w")
OUT_LIST.write ("genome\tlocus_tag\n")

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)

# Run HMM for a single marker
print("Searching for " + prot + " marker in genome set...")
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    outname= OUTPUT + "/out/"+dir+"/"+name + ".out"
    cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
    subprocess.call(cmd, stdout=FNULL)

# Parse HMM file 
print("Parsing all results...")
result_dir = os.walk(OUTPUT + "/out/"+dir)
for path, dirs, files in result_dir:
    for file in files:
        genome = file.replace(".out", "").strip().splitlines()[0]
        result = OUTPUT + "/out/"+dir+"/"+file
        output= OUTPUT + "/results/"+dir+".faa"
        genome_file=GENOMEDIR+genome+".faa"
        with open(output, "a") as outf:
            with open(genome_file, "r") as input_fasta:
                with open(result, "r") as input:
                    for qresult in SearchIO.parse(input, "hmmer3-tab"):
                        hits = qresult.hits
                        num_hits = len(hits)
                        if num_hits>0:
                            for i in range(0,1):
                                hit_id=hits[i].id
                            for record in SeqIO.parse(input_fasta, "fasta"):
                                if record.id in hit_id:
                                    outf.write(">"+genome+"\n"+str(record.seq)+"\n")
                                    OUT_LIST.write('%s\t%s\n' % (genome, record.id))
OUT_LIST.close()
outf.close()

# Align hits 
print("Aligning hits...")
prots = OUTPUT + "/results/*.faa"
fastas = glob.glob(prots)
for fasta in fastas:
    outname = os.path.basename(fasta).replace(".faa", "").strip().splitlines()[0]
    output= OUTPUT + "/results/"+outname+".aln"
    musc_cmd = ["muscle","-quiet","-in",fasta,"-out",output]
    subprocess.call(musc_cmd)

# Make tree 
if PHYTOOL == 'fastree':
    print("Calculating tree using FastTree...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    alignment_file = OUTPUT + "/results/"+marker_name+".aln"
    output_tree= OUTPUT + "/results/"+marker_name+".tre"
    tree_cmd = ["FastTree","-quiet","-out",output_tree,alignment_file]
    subprocess.call(tree_cmd)
elif PHYTOOL == "raxml":
    print("Calculating tree with RaxML... be patient...")
    marker_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
    outname= marker_name+"raxml"
    fileIn= OUTPUT + "/results/"+marker_name+".aln"
    raxCmd = "raxmlHPC-PTHREADS -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s "+fileIn+" -T "+THREADS+" -n "+outname
    os.system(raxCmd)

# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')