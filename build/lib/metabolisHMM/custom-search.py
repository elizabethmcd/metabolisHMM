#! /usr/bin/env python3 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# search-custom-markers to get statistics and visualization of a set of HMM markers
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import os, sys, glob, argparse, subprocess
import pandas as pd 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns
from distutils.spawn import find_executable
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO, SeqIO

# Arguments and Setup
parser = argparse.ArgumentParser(description = "Search custom directory of HMMs")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
kofam = parser.add_argument_group("kofam database arguments")
plotting = parser.add_argument_group("plotting arguments")

# required
required.add_argument('--input', metavar='GENOMEDIR', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help="Directory to store results and intermediate files")
required.add_argument('--markers_dir', metavar='MARKERDIR', help="Directory where custom markers are held")
required.add_argument('--markers_list', metavar='MARKERLIST', help="Ordered list of markers to run custom search on")
required.add_argument('--metadata', metavar='METADATA', help='Metadata file with taxonomical classifications or groups associated with genome file names')

# optional
optional.add_argument('--summary', metavar='OUTFILE', default='custom-markers-results.csv', help='Output statistics of custom marker searches in .csv format')
optional.add_argument('--heatmap', metavar='HEATOUT', default='custom-markers-results-heatmap.pdf', help="Summary heatmap of metabolic markers in PDF format. If you provide a custom name, it must end in .pdf" )

# kofam options
kofam.add_argument('--kofam', metavar='KOFAM', default = 'OFF', help="Use KofamKOALA HMM distributions. By default is OFF, options = ON,OFF")
kofam.add_argument("--ko_list", metavar='KOLIST', help="Point to location of the KofamKoala ko_list file if using the KofamKOALA KEGG HMMs")

# plotting options
plotting.add_argument('--aggregate', metavar='AGG', default='OFF', help="Aggregate metadata names by group = ON, visualize each genome individually = OFF" )
plotting.add_argument("--plotting", metavar='PLOTTING', default='ON', help="Option to turn plotting on or off depending on if you want just the raw statistics to perform your own plotting, or keep the default plotting settings. By default option is ON. Options= ON,OFF.")
plotting.add_argument("--ordering", metavar='ORDER', default='OFF', help="Turn ON or OFF to control custom ordering of rows in heatmap. By default is OFF to list in alphabetical ordering.")
plotting.add_argument("--group_list", metavar='GROUPS', help="Provide .txt file of ordering of groups if turn the aggregating option on to list the rows in custom order instead of the default alphabetical ordering.")

# check for required installed dependencies
# prodigal
if find_executable('prodigal') is not None:
    pass
else:
    print('You do not have prodigal installed in your path. Please fix this.')
    sys.exit()
# hmmsearch
if find_executable('hmmsearch') is not None:
    pass
else:
    print('You do not have the hmmsearch executable from HMMER in your path. Please fix this.')
    sys.exit()

# if no arguments given, print help message
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

# version to print
VERSION = '2.21'

# Beginning message
print('')
print('#############################################')
print('metabolisHMM v' + VERSION)


# arguments setup
args = parser.parse_args()
GENOMEDIR = args.input
GENOMEFILES = GENOMEDIR + "/**"
OUTPUT = args.output
OUTFILE = args.summary
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
out_genomes = OUTPUT + "/genomes"
KOFAM = args.kofam
KO_LIST = args.ko_list
PLOTTING = args.plotting
ORDER = args.ordering
GROUPS = args.group_list
MARKER_DIR = args.markers_dir
MARKERS_FILE = args.markers_list

# check if directory exists
if os.path.isdir(OUTPUT) == True:
    print("Directory "+ OUTPUT +" already exists! Please create different directory or remove the existing one.")
    sys.exit()

# check if directory of custom markers exists
if os.path.isdir(MARKER_DIR) == True:
    pass
    if os.path.exists(MARKERS_FILE) == False:
        print("You have not provided a list of the custom markers you would like run, and the order to put them in the heatmap. Please fix this.")
        sys.exit()
else:
    print("Directory for running custom HMM search was not provided. Please fix this.")
    sys.exit()

# if kofam argument turned on, check if the corresponding ko_list metadata is provided for getting TC scores from
if KOFAM == 'ON':
    if os.path.exists("ko_list") == False:
        print("     You have opted to use a marker from the KofamKOALA distribution, but have not supplied the path to the ko_list.")
        sys.exit()

# check contents of ko_list to make sure it's the correct file
if KOFAM == 'ON':
    with open(KO_LIST) as ko:
        first_line = ko.readlines()[0]
        elem = 'knum'
        if elem not in first_line:
            print("     This does not look like the ko_list file from the KofamKOALA distribution.")
            sys.exit()

# setup directories
os.makedirs(out_intm)
os.makedirs(out_results)
os.makedirs(out_genomes)
genomes = glob.glob(GENOMEFILES)

# point to directory where markers are, and only run on specific markers (in case want run on a select subset in a larger directory, such as the KofamKOALA database)
with open(MARKERS_FILE) as f:
    markers = f.read().splitlines()
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

print("Running HMM searches using custom marker set...")
for genome in reformatted_genomes: 
    name=os.path.basename(genome).replace(".reformatted.faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir(OUTPUT + "/out/"+dir)
    for marker in markers:
        infile=os.path.join(MARKER_DIR,marker)
        prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        outname= OUTPUT + "/out/"+dir+"/"+name + "-" + prot + ".out"
        if KOFAM == 'ON':
            kofam_name = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
            scores = pd.read_csv(KO_LIST, delimiter="\t")
            scoredic = scores.set_index('knum')['threshold'].to_dict()
            tc = scoredic.get(kofam_name)
            cmd = ["hmmsearch","-T"+tc,"--tblout="+outname,infile,genome]
            subprocess.call(cmd, stdout=FNULL)
        else: 
            cmd = ["hmmsearch","--cut_tc","--tblout="+outname, infile, genome]
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
df.fillna(0, inplace=True)

# Reformat dataframe in order of marker function, find markers/genomes with 0 hits and place 0's
# rows
all_rows = []
absent_rows = []
existing_rows = df.index
# columns
all_cols=[]
absent_cols=[]
existing_markers = df.columns
# check for missing rows
for genome in reformatted_genomes:
    name=os.path.basename(genome).replace(".reformatted.faa","")
    if name not in existing_rows:
        all_rows.append(name)
for row in existing_rows:
    all_rows.append(row)
df_rows = df.reindex(index=all_rows)
# check for missing cols
for marker in markers:
    prot=os.path.basename(marker).replace(".hmm","")
    if prot not in existing_markers:
        all_cols.append(prot)
for col in existing_markers:
    all_cols.append(col)
df_all = df_rows.reindex(columns=all_cols)
df_all.fillna(0, inplace=True)
# reorder columns
prots = []
with open(MARKERS_FILE) as f:
    markers = f.read().splitlines()
    for marker in markers:
        prot = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        prots.append(prot)
ordered_cols = prots
df_final = df_all[prots]
out_stats = OUTPUT + "/results/" + OUTFILE
df_final.to_csv(out_stats)

if PLOTTING == 'ON':
    print("Plotting results...")
    METADATA = args.metadata
    FIGOUT = args.heatmap
    FIGPATH = OUTPUT + "/results/" + FIGOUT
    AGG = args.aggregate
    # cleaned stats file
    clean = pd.read_csv(out_stats)
    clean.columns.values[0] = "genome"
    clean.set_index("genome", inplace=True)
    clean[clean>1] = 1
    out_clean = OUTPUT + "/results/cleaned-matrix.csv"
    clean.to_csv(out_clean)
    if AGG == 'OFF':
        if ORDER == 'OFF':
            sns.set(font_scale=2)
            fig, ax = plt.subplots(1, 1, figsize = (25, 15), dpi=300)
            xticks = clean.columns
            plot=sns.heatmap(clean, cmap="viridis",xticklabels=xticks, square=True, linewidths=1, linecolor='black', cbar=True, cbar_kws={"shrink": .50})
            plot.xaxis.tick_top()
            plot.xaxis.set_label_position('top')
            ax.set_ylabel('')
            plot.set_xticklabels(plot.get_xticklabels(), rotation=80)
            plot.figure.savefig(FIGPATH)
        elif ORDER == 'ON':
            if os.path.exists(GROUPS) == False:
                print("You have opted to list rows in custom order, but have not provided a groups list.")
                sys.exit()
            else: 
                rows = []
                with open(GROUPS) as f:
                    files = f.read().splitlines()
                    for file in files:
                        file = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
                        rows.append(file)
                clean.reindex(axis="index", level=0, labels=rows)
                sns.set(font_scale=2)
                fig, ax = plt.subplots(1, 1, figsize = (25, 15), dpi=300)
                xticks = clean.columns
                plot=sns.heatmap(clean, cmap="viridis",xticklabels=xticks, square=True, linewidths=1, linecolor='black', cbar=True, cbar_kws={"shrink": .50})
                plot.xaxis.tick_top()
                plot.xaxis.set_label_position('top')
                ax.set_ylabel('')
                plot.set_xticklabels(plot.get_xticklabels(), rotation=80)
                plot.figure.savefig(FIGPATH)
    elif AGG == 'ON':
        if ORDER == 'OFF':
            # merge with metadata
            metadata = pd.read_csv(METADATA, names=['genome', 'group'])
            merged = metadata.merge(clean, left_on="genome", right_on="genome")
            drop = merged.drop('genome', 1)
            agg = drop.groupby(['group'], as_index=False).mean()
            agg.set_index("group", inplace=True)
            sns.set(font_scale=2)
            fig, ax = plt.subplots(1, 1, figsize = (25, 15), dpi=300)
            xticks = agg.columns
            plot=sns.heatmap(agg, cmap="viridis",xticklabels=xticks, square=True, linewidths=1, linecolor='black', cbar=True, cbar_kws={"shrink": .50})
            plot.xaxis.tick_top()
            plot.xaxis.set_label_position('top')
            ax.set_ylabel('')
            plot.set_xticklabels(plot.get_xticklabels(), rotation=80)
            plot.figure.savefig(FIGPATH)
        elif ORDER == 'ON': 
            if os.path.exists(GROUPS) == False:
                print("You have opted to list rows in custom order, but have not provided a groups list.")
                sys.exit()
            else: 
                metadata = pd.read_csv(METADATA, names=['genome', 'group'])
                merged = metadata.merge(clean, left_on="genome", right_on="genome")
                drop = merged.drop('genome', 1)
                agg = drop.groupby(['group'], as_index=False).mean()
                agg.set_index("group", inplace=True)
                rows = []
                with open(GROUPS) as f:
                    files = f.read().splitlines()
                    for file in files:
                        file = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
                        rows.append(file)
                agg.reindex(axis="index", level=0, labels=rows)
                sns.set(font_scale=2)
                fig, ax = plt.subplots(1, 1, figsize = (25, 15), dpi=300)
                xticks = agg.columns
                plot=sns.heatmap(agg, cmap="viridis",xticklabels=xticks, square=True, linewidths=1, linecolor='black', cbar=True, cbar_kws={"shrink": .50})
                plot.xaxis.tick_top()
                plot.xaxis.set_label_position('top')
                ax.set_ylabel('')
                plot.set_xticklabels(plot.get_xticklabels(), rotation=80)
                plot.figure.savefig(FIGPATH)
        
elif PLOTTING == 'OFF':
    pass

# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')
