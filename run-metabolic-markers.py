#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
import pandas as pd 
from Bio import SearchIO

# Setup
genomes=glob.glob("genomes/*.faa")
markers=glob.glob("metabolic_markers/*.hmm")
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

# # Run HMMs
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    dir=name
    os.mkdir("out/"+dir)
    for marker in markers:
        prot=os.path.basename(marker).replace(".hmm", "").strip().splitlines()[0]
        outname= "out/"+dir+"/"+name + "-" + prot + ".out"
        cmd = ["hmmsearch","--cut_tc","--tblout="+outname, marker, genome]
        subprocess.call(cmd, stdout=FNULL)
        print("Running HMMsearch on " + name + " and " + prot + " marker")

# Parse HMM file to results matrix/dataframe
print("Parsing all results...")
all_dicts={}
result_dirs = os.walk("out/")
for path, dirs, files in result_dirs:
    genome_dict={}
    for file in files:
        genome = file.split("-")[0]
        prot = file.replace(".out", "").split("-")[1]
        result = "out/"+genome+"/"+file
        with open(result, "rU") as input:
            for qresult in SearchIO.parse(input, "hmmer3-tab"):
                hits = qresult.hits
                num_hits = len(hits)
                genome_dict[prot] = num_hits
                all_dicts[os.path.basename(file).split("-")[0]]=genome_dict
df=pd.DataFrame.from_dict(all_dicts, orient="index", dtype=None)
df.to_csv("results/metabolic-results.tsv", sep="\t")

# Reformat dataframe in order of marker function, find markers in none of the genomes and input NaN for all
absent_cols=[]
existing_markers = df.columns
all_markers=glob.glob("metabolic_markers/*.hmm")
for marker in all_markers:
    prot=os.path.basename(marker).replace(".hmm", "")
    if prot not in existing_markers:
        absent=prot
        absent_cols.append(absent)
print(absent_cols)

reordered_cols = ['fccB', 'sulfide_quinone_oxidoreductase_sqr', 'sulfur_dioxygenase_sdo', 'aprA', 'sat', 'cysC', 'cysN', 'thiosulfate_reductase_phsA', 'FeFeHydrogenase', 'FeFeHydrogenase_2', 'Hydrogenase_Group_1', 'Hydrogenase_Group_2a', 'Hydrogenase_Group_2b', 'Hydrogenase_Group_3b', 'Hydrogenase_Group_3c', 'Hydrogenase_Group_3d', 'Hydrogenase_Group_4', 'pmoA', 'pmoB', 'pmoC', 'mmoB', 'mmoD', 'mcrA', 'mcrB', 'mcrC', 'nifA_Mo', 'nifB_Mo', 'nifH', 'nitrite_oxidoreductase_nxrA', 'nitrite_oxidoreductase_nxrB', 'napA', 'napB', 'narG', 'narH', 'nrfH', 'nrfA', 'nirB', 'nirD', 'nirK', 'nitrite_reductase_nirS', 'nitric_oxide_reductase_norB', 'nitric_oxide_reductase_norC', 'nosD', 'nosZ', 'hydrazine_oxidase_hzoA', 'hydrazine_synthase_hzsA', 'octR', 'coxA', 'coxB', 'ccoN', 'ccoO', 'ccoP', 'cyoA', 'cyoD', 'cyoE', 'cydA', 'cydB', 'qoxA', 'ndma_methanol_dehydrogenase', 'madA', 'madB', 'fdh_thiol_id', 'sfh', 'sgdh', 'smdh', 'fae', 'fmtf', 'mtmc', 'fdhA', 'fdhB', 'fdhC', 'carbon_monoxide_dehydrogenase_coxL', 'carbon_monoxide_dehydrogenase_coxM', 'carbon_monoxide_dehydrogenase_coxS', 'rubisco_form_I', 'rubisco_form_II_III', 'rubisco_form_III', 'rubisco_form_IV', 'propionyl_coA_synthase', 'malonyl_coA_reductase', '4_hydroxybutyryl_coA_dehydratase', '4_hydroxybutyryl_coA_synthetase', 'codhC', 'codhD', 'codh_catalytic', 'acetate_citrate_lyase_aclA', 'acetate_citrate_lyase_aclB', 'ureA', 'ureB', 'ureC', 'hdh', 'rdh', 'cld', 'ars_ox', 'ars_ox_2', 'arsC', 'ars_thioredoxin', 'ars_glutaredoxin', 'ygfK', 'ygfM', 'sel_mo', 'nthA', 'nthB', 'mtrA', 'mtrC', 'dsrA', 'dsrB', 'dsrC', 'dsrD', 'dsrE', 'dsrF', 'dsrH', 'dsrJ', 'dsrK', 'dsrM', 'dsrO', 'dsrP', 'dsrR', 'dsrS', 'dsrT', 'soxA', 'soxB', 'soxD', 'soxX', 'soxY', 'soxYZ-like', 'soxZ']


