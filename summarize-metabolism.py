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
required.add_argument('--output', metavar='OUTPUT', help="Directory to sotre results and intermediate files")
optional.add_argument('--summary', metavar='OUTFILE', default="metabolic-summary-results.csv", help="Summary file of metabolic marker statistics in CSV format")
optional.add_argument('--heatmap', metavar='HEATOUT', default='metabolic-summary-results-heatmap.pdf', help="Summary heatmap of metabolic markers" )
optional.add_argument('--metadata', metavar='METADATA', help='Metadata file with taxonomical classifications or groups associated with genome file names')
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

# arguments and directory setup
args = parser.parse_args()
GENOMEDIR = args.input
OUTPUT = args.output
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
OUTFILE = args.summary
genomes=glob.glob(os.path.join(GENOMEDIR, '*.faa'))
markers=glob.glob("metabolic_markers/*.hmm")
os.makedirs(out_intm)
os.makedirs(out_results)
FNULL = open(os.devnull, 'w')

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)

# Run HMMs
print("Running metabolic marker HMMsearches...")
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
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



# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')
