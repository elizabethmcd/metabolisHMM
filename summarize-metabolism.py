#! /usr/bin/env python 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# summarize-metabolism : get statistics and visualizations of common metabolic markers spanning multiple biogeochemical cycles
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import os, sys, glob, subprocess, argparse
import pandas as pd 
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

# Arguments and Setup 
parser = argparse.ArgumentParser(description = "Create a summary of metabolic markers spanning biogeochemical cycles")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
required.add_argument('--input', metavar='GENOMEDIR', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help="Directory to store results and intermediate files")
optional.add_argument('--summary', metavar='OUTFILE', default="metabolic-summary-results.csv", help="Summary file of metabolic marker statistics in CSV format")
optional.add_argument('--heatmap', metavar='HEATOUT', default='metabolic-summary-results-heatmap.pdf', help="Summary heatmap of metabolic markers in PDF format. If you provide a custom name, it must end in .pdf" )
required.add_argument('--metadata', metavar='METADATA', help='Metadata file with taxonomical classifications or groups associated with genome file names. Needs to be a .csv file with the genome/filename (everything before .fna) and the group name.')
optional.add_argument('--aggregate', metavar='AGG', default='OFF', help="Aggregate metadata names by group = ON, visualize each genome individually = OFF" )

# if no arguments given, print help message
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# version to print
def version():
    versionFile = open('VERSION')
    return versionFile.read()
VERSION = version()

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)

# arguments setup
args = parser.parse_args()
GENOMEDIR = args.input
GENOMEFILES = GENOMEDIR + "/**"
OUTPUT = args.output
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
out_genomes = OUTPUT + "/genomes"
OUTFILE = args.summary
markers=glob.glob("metabolic_markers/*.hmm")
genomes=glob.glob(GENOMEFILES)

# check if directory exists
if os.path.isdir(OUTPUT) == True:
    print("Directory "+ OUTPUT +" already exists! Please create different directory or remove the existing one.")
    sys.exit()

# make directories
os.makedirs(out_intm)
os.makedirs(out_results)
os.makedirs(out_genomes)
# this also fixed with subprocess DEVNULL
FNULL = open(os.devnull, 'w')

# if .fna predict CDS and reformat header names because prodigal makes them stupid
# if .faa reformat the headers just in case contains weirdness
# if the user didn't provide the right files tell them
n = 0
print("Reformatting fasta files...")
for genome in genomes:
    if genome.endswith('.fna'):
        name = os.path.basename(genome).replace(".fna", "").strip().splitlines()[0]
        out_prot = OUTPUT + "/genomes/" + name + ".faa"
        out_gbk = OUTPUT + "/genomes/" + name + ".gbk"
        out_reformatted = OUTPUT + "/genomes/" + name + ".reformatted.faa"
        prodigal_cmd = "prodigal -q -i "+genome+" -a "+out_prot +" -o "+out_gbk
        os.system(prodigal_cmd)
        for seq_record in SeqIO.parse(out_prot, "fasta"):
            n = n + 1
            a = str(n).zfill(5)
            with open(out_reformatted, "a") as outre:
                outre.write(">" + name + "_" + str(a) + "\n")
                outre.write(str(seq_record.seq) + "\n")
    elif genome.endswith('.faa'):
        name = os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
        out_reformatted = OUTPUT + "/genomes/" + name + ".reformatted.faa"
        for seq_record in SeqIO.parse(genome, "fasta"):
            n = n + 1
            a = str(n).zfill(5)
            with open(out_reformatted, "a") as outre:
                outre.write(">" + name + "_" + str(a) + "\n")
                outre.write(str(seq_record.seq) + "\n")
    else:
        print("These do not look like fasta files that end in .fna or .faa. Please check your genome files.")
        sys.exit()
reformatted_path = OUTPUT + "/genomes/" + "*.reformatted.faa"
reformatted_genomes = glob.glob(reformatted_path)

# Run HMMs
print("Running metabolic marker HMMsearches...")
for genome in reformatted_genomes: 
    name=os.path.basename(genome).replace(".reformatted.faa", "").strip().splitlines()[0]
    dir=name
    os.makedirs(OUTPUT + "/out/"+dir)
    for marker in markers:
        prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        outname= OUTPUT + "/out/"+dir+"/"+name + "-" + prot + ".out"
        cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
        subprocess.call(cmd, stdout=FNULL)

# Parse HMM file to results matrix/dataframe
print("Parsing all results...")
all_dicts={}
result_dirs = os.walk(OUTPUT + "/out/")
for path, dirs, files in result_dirs:
    genome_dict={}
    for file in files:
        genome = file.split("-")[0]
        prot = file.replace(".out", "").split("-")[1]
        result = OUTPUT + "/out/"+genome+"/"+file
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

# open list in data/metabolic-markers.txt to reorder by function
with open("data/marker-list.txt") as f:
    reordered_cols = f.read().splitlines()
df_final=df_all[reordered_cols]
out_stats = OUTPUT + "/results/" + OUTFILE
df_final.to_csv(out_stats)

print("Plotting results...")
METADATA = args.metadata
FIGOUT = args.heatmap
AGG = args.aggregate
plot_cmd = ['Rscript','make-heatmap.R', out_stats, METADATA, AGG, FIGOUT]
subprocess.call(plot_cmd, stdout=FNULL)

# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')
