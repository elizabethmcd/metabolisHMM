#! /usr/bin/env python3 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# summarize-metabolism : get statistics and visualizations of common metabolic markers spanning multiple biogeochemical cycles
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import os, sys, glob, subprocess, argparse
import pandas as pd 
import numpy as np
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
parser = argparse.ArgumentParser(description = "Create a summary of metabolic markers spanning biogeochemical cycles")
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
plotting = parser.add_argument_group("plotting arguments")

# required
required.add_argument('--input', metavar='GENOMEDIR', help='Directory where genomes to be screened are held')
required.add_argument('--output', metavar='OUTPUT', help="Directory to store results and intermediate files")
required.add_argument('--metadata', metavar='METADATA', help='Metadata file with taxonomical classifications or groups associated with genome file names. Needs to be a .csv file with the genome/filename (everything before .fna) and the group name.')

# optional
optional.add_argument('--summary', metavar='OUTFILE', default="raw-metabolic-summary-results.csv", help="Summary file of metabolic marker statistics in CSV format")
optional.add_argument('--heatmap', metavar='HEATOUT', default='metabolic-summary-results-heatmap.pdf', help="Summary heatmap of metabolic markers in PDF format. If you provide a custom name, it must end in .pdf" )

# plotting
optional.add_argument('--aggregate', metavar='AGG', default='OFF', help="Aggregate metadata names by group = ON, visualize each genome individually = OFF" )
optional.add_argument("--plotting", metavar='PLOTTING', default='ON', help="Option to turn plotting on or off depending on if you want just the raw statistics to perform your own plotting, or keep the default plotting settings. Options= ON,OFF.")
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
out_intm = OUTPUT + "/out"
out_results = OUTPUT + "/results"
out_genomes = OUTPUT + "/genomes"
OUTFILE = args.summary
markers=glob.glob("curated_markers/metabolic_markers/*.hmm")
genomes=glob.glob(GENOMEFILES)
PLOTTING = args.plotting

# check if directory exists
if os.path.isdir(OUTPUT) == True:
    print("Directory "+ OUTPUT +" already exists! Please create different directory or remove the existing one.")
    sys.exit()

# check for curated metabolic markers direcotry
if os.path.isdir("curated_markers/metabolic_markers/") == False:
    print("     The directory of curated metabolic markers could not be found."+"\n"+"     Please either download the markers from https://github.com/elizabethmcd/metabolisHMM/releases/download/v2.0/metabolisHMM_v2.0_markers.tgz and decompress the tarball, or move the directory to where you are running the workflow from.")
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
print("Screening curated metabolic markers...")
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

# Reformat dataframe in order of marker function, find markers/genomes with no results and place 0's
# rows
all_rows = []
absent_rows = []
existing_rows = df.index
# columns
all_cols=[]
absent_cols=[]
existing_markers = df.columns
markers=glob.glob("curated_markers/metabolic_markers/*.hmm")
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

# open list in data/metabolic-markers.txt to reorder by function
reordered_cols = ['ndma_methanol_dehydrogenase_TIGR04266','madA_TIGR02659','madB_TIGR02658','fdh_thiol_id_TIGR02819','sfh_TIGR02821','sgdh_TIGR02818','smdh_TIGR03451','fae_TIGR03126','fmtf_TIGR03119','mtmc_TIGR03120','fdhA_TIGR01591','fdhB_TIGR01582','fdhC_TIGR01583','carbon_monoxide_dehydrogenase_coxL_TIGR02416','carbon_monoxide_dehydrogenase_coxM','carbon_monoxide_dehydrogenase_coxS','rubisco_form_I','rubisco_form_II','rubisco_form_III','rubisco_form_II_III','rubisco_form_IV','codhC_TIGR00316','codhD_TIGR00381','codh_catalytic_TIGR01702','acetate_citrate_lyase_aclA','acetate_citrate_lyase_aclB','nifD_TIGR01282','nifH_TIGR01287','nifK_TIGR01286','nitrite_oxidoreductase_nxrA','nitrite_oxidoreductase_nxrB','napA_TIGR01706','napB_PF03892','narG_TIGR01580','narH_TIGR01660','nrfA_TIGR03152','nrfH_TIGR03153','nirB_TIGR02374','nirD_TIGR02378','nirK_TIGR02376','nitrite_reductase_nirS','nitric_oxide_reductase_norB','nitric_oxide_reductase_norC','nosD_TIGR04247','nosZ_TIGR04246','hydrazine_oxidase_hzoA','hydrazine_synthase_hzsA','fccB_PF09242','sulfide_quinone_oxidoreductase_sqr','sulfur_dioxygenase_sdo','aprA_TIGR02061','sat_TIGR00339','dsrA_TIGR02064','dsrB_TIGR02066','dsrD_PF08679','thiosulfate_reductase_phsA','soxB_TIGR04486','soxC_TIGR04555','soxY_TIGR04488','coxA_TIGR02891','coxB_TIGR02866','ccoN_TIGR00780','ccoO_TIGR00781','ccoP_TIGR00782','cyoA_TIGR01433','cyoD_TIGR02847','cyoE_TIGR01473','cydA_PF01654','cydB_TIGR00203','qoxA_TIGR01432','FeFeHydrogenase_TIGR02512','FeFeHydrogenase_TIGR04105','Hydrogenase_Group_1','Hydrogenase_Group_2a','Hydrogenase_Group_2b','Hydrogenase_Group_3a','Hydrogenase_Group_3b','Hydrogenase_Group_3c','Hydrogenase_Group_3d','Hydrogenase_Group_4']

df_final=df_all[reordered_cols]
out_stats = OUTPUT + "/results/" + OUTFILE
df_final.to_csv(out_stats)

if PLOTTING == 'ON':
    METADATA = args.metadata
    FIGOUT = args.heatmap
    FIGPATH = OUTPUT + "/results/" + FIGOUT
    AGG = args.aggregate
    ORDER = args.ordering
    GROUPS = args.group_list
    clean = pd.read_csv(out_stats)
    clean = pd.read_csv(out_stats)
    clean.columns.values[0] = "genome"
    clean.set_index("genome", inplace=True)
    clean[clean>1] = 1
    out_clean = OUTPUT + "/results/cleaned-matrix.csv"
    clean.to_csv(out_clean)
    clean.reset_index(inplace=True)
    if AGG == 'OFF':
        if ORDER == 'OFF':
            sns.set(font_scale=2)
            carbon = clean.iloc[:,0:27].copy()
            carbon.columns = ["genome","MtOH dehy", "madA", "madB", "fdh", "sfh", "sgdh", "smdh", "fae", "fmtF", "mtmc", "fdhA", "fdhB", "fdhC", "coxL", "coxM", "coxS", "rubisco I", "rubisco II", "rubisco III", "rubisco II/III", "rubisco IV", "codhC", "codhD", "codh cat.", "aclA", "aclB"]
            carbon.set_index("genome", inplace=True)
            nitrogen = clean.iloc[:, np.r_[0, 27:48]].copy()
            nitrogen.columns= ["genome","nifD", "nifH", "nifK", "nxrA", "nxrB", "napA", "napB", "narG", "narH", "nrfA", "nrfH", "nirB", "nirD", "nirK", "nirS", "norB", "norC", "nosD", "nosZ", "hzoA", "hzsA"]
            nitrogen.set_index("genome", inplace=True)
            sulfur = clean.iloc[:, np.r_[0,48:60]].copy()
            sulfur.columns= ["genome","fccB", "sqr", "sdo", "aprA", "sat", "dsrA", "dsrB", "dsrD", "phsA", "soxB", "soxC", "soxY"]
            sulfur.set_index("genome", inplace=True)
            oxygen = clean.iloc[:, np.r_[0,60:71]].copy()
            oxygen.columns = ["genome","coxA", "coxB", "ccoN", "ccoO", "ccoP", "cyoA", "cyoD", "cyoE", "cydA", "cydB", "qoxA"]
            oxygen.set_index("genome", inplace=True)
            hydrogen = clean.iloc[:, np.r_[0,71:81]].copy()
            hydrogen.columns = ["genome","FeFe I", "FeFe II", "Group I", "Group IIA", "Group IIB", "Group IIIA", "Group IIIB", "Group IIIC", "Group IIID", "Group IV"]
            hydrogen.set_index("genome", inplace=True)
            fig, (ax,ax2,ax3,ax4,ax5) = plt.subplots(ncols=5, figsize = (85, 20))
            fig.subplots_adjust(wspace=0.01)
            sns.heatmap(carbon, cmap="Blues", ax=ax, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(nitrogen, cmap="Greens", ax=ax2, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(sulfur, cmap="Purples", ax=ax3, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(oxygen, cmap="Greys", ax=ax4, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(hydrogen, cmap="Oranges", ax=ax5, cbar=False, linewidths=1, linecolor='black')
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=80)
            ax2.set_ylabel('')
            ax2.xaxis.tick_top()
            ax2.xaxis.set_label_position('top')
            ax2.set_yticks([])
            ax2.tick_params(rotation=80)
            ax3.set_ylabel('')
            ax3.xaxis.tick_top()
            ax3.xaxis.set_label_position('top')
            ax3.set_yticks([])
            ax3.tick_params(rotation=80)
            ax4.set_ylabel('')
            ax4.xaxis.tick_top()
            ax4.xaxis.set_label_position('top')
            ax4.set_yticks([])
            ax4.tick_params(rotation=80)
            ax5.xaxis.tick_top()
            ax5.xaxis.set_label_position('top')
            ax5.set_ylabel('')
            ax5.set_yticks([])
            ax5.tick_params(rotation=80)
            plt.savefig(FIGPATH)
        elif ORDER == 'ON':
            if os.path.exists(GROUPS) == False:
                print("You have opted to list rows in custom order, but have not provided a groups list.")
                sys.exit()
            else: 
                rows = []
                with open(METADATA) as f:
                    files = f.read().splitlines()
                    for file in files:
                        file = os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
                        rows.append(file)
                clean.reindex(axis="index", level=0, labels=rows)
                clean.reset_index()
                sns.set(font_scale=2)
                carbon = clean.iloc[:,0:27].copy()
                carbon.columns = ["group","MtOH dehy", "madA", "madB", "fdh", "sfh", "sgdh", "smdh", "fae", "fmtF", "mtmc", "fdhA", "fdhB", "fdhC", "coxL", "coxM", "coxS", "rubisco I", "rubisco II", "rubisco III", "rubisco II/III", "rubisco IV", "codhC", "codhD", "codh cat.", "aclA", "aclB"]
                carbon.set_index("group", inplace=True)
                nitrogen = clean.iloc[:, np.r_[0, 27:48]].copy()
                nitrogen.columns= ["group","nifD", "nifH", "nifK", "nxrA", "nxrB", "napA", "napB", "narG", "narH", "nrfA", "nrfH", "nirB", "nirD", "nirK", "nirS", "norB", "norC", "nosD", "nosZ", "hzoA", "hzsA"]
                nitrogen.set_index("group", inplace=True)
                sulfur = clean.iloc[:, np.r_[0,48:60]].copy()
                sulfur.columns= ["group","fccB", "sqr", "sdo", "aprA", "sat", "dsrA", "dsrB", "dsrD", "phsA", "soxB", "soxC", "soxY"]
                sulfur.set_index("group", inplace=True)
                oxygen = clean.iloc[:, np.r_[0,60:71]].copy()
                oxygen.columns = ["group","coxA", "coxB", "ccoN", "ccoO", "ccoP", "cyoA", "cyoD", "cyoE", "cydA", "cydB", "qoxA"]
                oxygen.set_index("group", inplace=True)
                hydrogen = clean.iloc[:, np.r_[0,71:81]].copy()
                hydrogen.columns = ["group","FeFe I", "FeFe II", "Group I", "Group IIA", "Group IIB", "Group IIIA", "Group IIIB", "Group IIIC", "Group IIID", "Group IV"]
                hydrogen.set_index("group", inplace=True)
                fig, (ax,ax2,ax3,ax4,ax5) = plt.subplots(ncols=5, figsize = (85, 20))
                fig.subplots_adjust(wspace=0.01)
                sns.heatmap(carbon, cmap="Blues", ax=ax, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(nitrogen, cmap="Greens", ax=ax2, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(sulfur, cmap="Purples", ax=ax3, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(oxygen, cmap="Greys", ax=ax4, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(hydrogen, cmap="Oranges", ax=ax5, cbar=False, linewidths=1, linecolor='black')
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
                ax.set_xticklabels(ax.get_xticklabels(), rotation=80)
                ax2.set_ylabel('')
                ax2.xaxis.tick_top()
                ax2.xaxis.set_label_position('top')
                ax2.set_yticks([])
                ax2.tick_params(rotation=80)
                ax3.set_ylabel('')
                ax3.xaxis.tick_top()
                ax3.xaxis.set_label_position('top')
                ax3.set_yticks([])
                ax3.tick_params(rotation=80)
                ax4.set_ylabel('')
                ax4.xaxis.tick_top()
                ax4.xaxis.set_label_position('top')
                ax4.set_yticks([])
                ax4.tick_params(rotation=80)
                ax5.xaxis.tick_top()
                ax5.xaxis.set_label_position('top')
                ax5.set_ylabel('')
                ax5.set_yticks([])
                ax5.tick_params(rotation=80)
                plt.savefig(FIGPATH)
    elif AGG == 'ON':
        if ORDER == 'OFF':
            metadata = pd.read_csv(METADATA, names=['genome', 'group'])
            merged = metadata.merge(clean, left_on="genome", right_on="genome")
            drop = merged.drop('genome', 1)
            agg = drop.groupby(['group'], as_index=False).mean()
            sns.set(font_scale=3)
            carbon = agg.iloc[:,0:27].copy()
            carbon.columns = ["group","MtOH dehy", "madA", "madB", "fdh", "sfh", "sgdh", "smdh", "fae", "fmtF", "mtmc", "fdhA", "fdhB", "fdhC", "coxL", "coxM", "coxS", "rubisco I", "rubisco II", "rubisco III", "rubisco II/III", "rubisco IV", "codhC", "codhD", "codh cat.", "aclA", "aclB"]
            carbon.set_index("group", inplace=True)
            nitrogen = agg.iloc[:, np.r_[0, 27:48]].copy()
            nitrogen.columns= ["group","nifD", "nifH", "nifK", "nxrA", "nxrB", "napA", "napB", "narG", "narH", "nrfA", "nrfH", "nirB", "nirD", "nirK", "nirS", "norB", "norC", "nosD", "nosZ", "hzoA", "hzsA"]
            nitrogen.set_index("group", inplace=True)
            sulfur = agg.iloc[:, np.r_[0,48:60]].copy()
            sulfur.columns= ["group","fccB", "sqr", "sdo", "aprA", "sat", "dsrA", "dsrB", "dsrD", "phsA", "soxB", "soxC", "soxY"]
            sulfur.set_index("group", inplace=True)
            oxygen = agg.iloc[:, np.r_[0,60:71]].copy()
            oxygen.columns = ["group","coxA", "coxB", "ccoN", "ccoO", "ccoP", "cyoA", "cyoD", "cyoE", "cydA", "cydB", "qoxA"]
            oxygen.set_index("group", inplace=True)
            hydrogen = agg.iloc[:, np.r_[0,71:81]].copy()
            hydrogen.columns = ["group","FeFe I", "FeFe II", "Group I", "Group IIA", "Group IIB", "Group IIIA", "Group IIIB", "Group IIIC", "Group IIID", "Group IV"]
            hydrogen.set_index("group", inplace=True)
            fig, (ax,ax2,ax3,ax4,ax5) = plt.subplots(ncols=5, figsize = (85, 20))
            fig.subplots_adjust(wspace=0.01)
            sns.heatmap(carbon, cmap="Blues", ax=ax, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(nitrogen, cmap="Greens", ax=ax2, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(sulfur, cmap="Purples", ax=ax3, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(oxygen, cmap="Greys", ax=ax4, cbar=False, linewidths=1, linecolor='black')
            sns.heatmap(hydrogen, cmap="Oranges", ax=ax5, cbar=False, linewidths=1, linecolor='black')
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=80)
            ax2.set_ylabel('')
            ax2.xaxis.tick_top()
            ax2.xaxis.set_label_position('top')
            ax2.set_yticks([])
            ax2.tick_params(rotation=80)
            ax3.set_ylabel('')
            ax3.xaxis.tick_top()
            ax3.xaxis.set_label_position('top')
            ax3.set_yticks([])
            ax3.tick_params(rotation=80)
            ax4.set_ylabel('')
            ax4.xaxis.tick_top()
            ax4.xaxis.set_label_position('top')
            ax4.set_yticks([])
            ax4.tick_params(rotation=80)
            ax5.xaxis.tick_top()
            ax5.xaxis.set_label_position('top')
            ax5.set_ylabel('')
            ax5.set_yticks([])
            ax5.tick_params(rotation=80)
            plt.savefig(FIGPATH)
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
                agg.reset_index()
                sns.set(font_scale=3)
                carbon = agg.iloc[:,0:27].copy()
                carbon.columns = ["group","MtOH dehy", "madA", "madB", "fdh", "sfh", "sgdh", "smdh", "fae", "fmtF", "mtmc", "fdhA", "fdhB", "fdhC", "coxL", "coxM", "coxS", "rubisco I", "rubisco II", "rubisco III", "rubisco II/III", "rubisco IV", "codhC", "codhD", "codh cat.", "aclA", "aclB"]
                carbon.set_index("group", inplace=True)
                nitrogen = agg.iloc[:, np.r_[0, 27:48]].copy()
                nitrogen.columns= ["group","nifD", "nifH", "nifK", "nxrA", "nxrB", "napA", "napB", "narG", "narH", "nrfA", "nrfH", "nirB", "nirD", "nirK", "nirS", "norB", "norC", "nosD", "nosZ", "hzoA", "hzsA"]
                nitrogen.set_index("group", inplace=True)
                sulfur = agg.iloc[:, np.r_[0,48:60]].copy()
                sulfur.columns= ["group","fccB", "sqr", "sdo", "aprA", "sat", "dsrA", "dsrB", "dsrD", "phsA", "soxB", "soxC", "soxY"]
                sulfur.set_index("group", inplace=True)
                oxygen = agg.iloc[:, np.r_[0,60:71]].copy()
                oxygen.columns = ["group","coxA", "coxB", "ccoN", "ccoO", "ccoP", "cyoA", "cyoD", "cyoE", "cydA", "cydB", "qoxA"]
                oxygen.set_index("group", inplace=True)
                hydrogen = agg.iloc[:, np.r_[0,71:81]].copy()
                hydrogen.columns = ["group","FeFe I", "FeFe II", "Group I", "Group IIA", "Group IIB", "Group IIIA", "Group IIIB", "Group IIIC", "Group IIID", "Group IV"]
                hydrogen.set_index("group", inplace=True)
                fig, (ax,ax2,ax3,ax4,ax5) = plt.subplots(ncols=5, figsize = (85, 20))
                fig.subplots_adjust(wspace=0.01)
                sns.heatmap(carbon, cmap="Blues", ax=ax, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(nitrogen, cmap="Greens", ax=ax2, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(sulfur, cmap="Purples", ax=ax3, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(oxygen, cmap="Greys", ax=ax4, cbar=False, linewidths=1, linecolor='black')
                sns.heatmap(hydrogen, cmap="Oranges", ax=ax5, cbar=False, linewidths=1, linecolor='black')
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top')
                ax.set_xticklabels(ax.get_xticklabels(), rotation=80)
                ax2.set_ylabel('')
                ax2.xaxis.tick_top()
                ax2.xaxis.set_label_position('top')
                ax2.set_yticks([])
                ax2.tick_params(rotation=80)
                ax3.set_ylabel('')
                ax3.xaxis.tick_top()
                ax3.xaxis.set_label_position('top')
                ax3.set_yticks([])
                ax3.tick_params(rotation=80)
                ax4.set_ylabel('')
                ax4.xaxis.tick_top()
                ax4.xaxis.set_label_position('top')
                ax4.set_yticks([])
                ax4.tick_params(rotation=80)
                ax5.xaxis.tick_top()
                ax5.xaxis.set_label_position('top')
                ax5.set_ylabel('')
                ax5.set_yticks([])
                ax5.tick_params(rotation=80)
                plt.savefig(FIGPATH)

elif PLOTTING == 'OFF':
    pass

# end message
print("Done! Find your results in "+ OUTPUT + "/results/")
print('#############################################')
