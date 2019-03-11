#! /usr/bin/env python 

import os, sys
import glob 
import subprocess 
import pandas as pd 
from Bio import SearchIO

# Usage 
parser = argparse.ArgumentParser(description = "Output summary of provided metabolic markers against input genomes")

# Setup
genomes=glob.glob("genomes/*.faa")
markers=glob.glob("metabolic_markers/*.hmm")
os.mkdir("out")
os.mkdir("results")
FNULL = open(os.devnull, 'w')

# Run HMMs
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

reordered_cols = ['fccB_PF09242', 'sulfide_quinone_oxidoreductase_sqr', 'sor_PF07682', 'sulfur_dioxygenase_sdo', 'aprA_TIGR02061', 'sat_TIGR00339', 'cysC_TIGR00455', 'cysN_TIGR02034', 'thiosulfate_reductase_phsA', 'asrA_TIGR02910', 'asrB_TIGR02911', 'asrC_TIGR02912', 'FeFeHydrogenase_TIGR02512', 'FeFeHydrogenase_TIGR04105', 'Hydrogenase_Group_1', 'Hydrogenase_Group_2a', 'Hydrogenase_Group_2b', 'Hydrogenase_Group_3a', 'Hydrogenase_Group_3b', 'Hydrogenase_Group_3c', 'Hydrogenase_Group_3d', 'Hydrogenase_Group_4', 'pmoA_TIGR03080', 'pmoB_TIGR03079', 'pmoC_TIGR03078', 'mmoB_PF02406', 'mmoD_TIGR04550', 'mcrA_TIGR03256', 'mcrB_TIGR03257', 'mcrC_TIGR03259', 'nifA_Mo_TIGR01282', 'nifB_Mo_TIGR01286', 'nifH_TIGR01287', 'nifA_Fe_TIGR01861', 'nifA_Va_TIGR01860', 'nifB_Fe_TIGR02931', 'nifB_V_TIGR02932', 'nifD_Fe_TIGR02929', 'nifD_V_TIGR02930', 'nitrite_oxidoreductase_nxrA', 'nitrite_oxidoreductase_nxrB', 'napA_TIGR01706', 'napB_PF03892', 'narG_TIGR01580', 'narH_TIGR01660', 'nrfH_TIGR03153', 'nrfA_PF02335', 'nrfA_TIGR03152', 'nrfD_TIGR03148', 'nirB_TIGR02374', 'nirD_TIGR02378', 'nirK_TIGR02376', 'nitrite_reductase_nirS', 'nitric_oxide_reductase_norB', 'nitric_oxide_reductase_norC', 'nosD_TIGR04247', 'nosZ_TIGR04246', 'hydrazine_oxidase_hzoA', 'hydrazine_synthase_hzsA', 'octR_TIGR04315', 'coxA_TIGR02891', 'coxB_TIGR02866', 'ccoN_TIGR00780', 'ccoO_TIGR00781', 'ccoP_TIGR00782', 'cyoA_TIGR01433', 'cyoD_TIGR02847', 'cyoE_TIGR01473', 'cydA_PF01654', 'cydB_TIGR00203', 'qoxA_TIGR01432', 'qoxB_TIGR02882', 'ndma_methanol_dehydrogenase_TIGR04266', 'methanol_dehydrogenase_pqq_xoxF_mxaF', 'madA_TIGR02659', 'madB_TIGR02658', 'fdh_thiol_id_TIGR02819', 'sfh_TIGR02821', 'sgdh_TIGR02818', 'smdh_TIGR03451', 'fae_TIGR03126', 'fmtf_TIGR03119', 'mtmc_TIGR03120', 'fdhA_TIGR01591', 'fdhB_TIGR01582', 'fdhC_TIGR01583', 'carbon_monoxide_dehydrogenase_coxL_TIGR02416', 'carbon_monoxide_dehydrogenase_coxM', 'carbon_monoxide_dehydrogenase_coxS', 'rubisco_form_I', 'rubisco_form_II', 'rubisco_form_III', 'rubisco_form_II_III', 'rubisco_form_IV', 'propionyl_coA_synthase', 'malonyl_coA_reductase', '4_hydroxybutyryl_coA_dehydrogenase', '4_hydroxybutyryl_coA_synthetase', 'codhC_TIGR00316', 'codhD_TIGR00381', 'codh_catalytic_TIGR01702', 'acetate_citrate_lyase_aclA', 'acetate_citrate_lyase_aclB', 'ureA_TIGR00193', 'ureB_TIGR00192', 'ureC_TIGR01792', 'hdh_TIGR01428', 'rdh_TIGR02486', 'cld_PF06778', 'ars_ox_TIGR02693', 'ars_ox_TIGR02694', 'arsC_TIGR00014', 'ars_thioredoxin_TIGR02691', 'ars_glutaredoxin_TIGR02689', 'ygfK_TIGR03315', 'ygfM_TIGR03312', 'sel_mo_TIGR03313', 'nthA_TIGR01323', 'nthB_TIGR03888', 'mtrA_TIGR03507', 'mtrC_TIGR03509', 'dsrA_TIGR02064', 'dsrB_TIGR02066', 'dsrD_PF08679','soxB_TIGR04486','soxC_TIGR04555','soxY_TIGR04488', 'pcrA_TIGR03479', 'pcrB_TIGR03478', 'hgcA']

df_final=df_all[reordered_cols]
df_final.to_csv("results/metabolic-marker-results.tsv", sep="\t")
df_final.to_csv("results/metabolic-marker-results.csv")
