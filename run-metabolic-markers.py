#! /usr/bin/env python 

import os, sys
import glob 
from Bio import SearchIO

# Setup
genomes=glob.glob("genomes/*.faa")
markers=glob.glob("metabolic_markers/*.hmm")
FNULL = open(os.devnull, 'w')

# Run HMMs
for genome in genomes: 
    name=os.path.basename(genome).replace(".faa", "").strip().splitlines()[0]
    for marker in markers:
        