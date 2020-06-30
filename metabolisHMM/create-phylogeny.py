#! /usr/bin/env python3 

###############################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# create-genome-phylogeny: to create a phylogeny based on ribosomal proteins
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
###############################################################

import glob, argparse, subprocess, os, sys, tempfile, re
from subprocess import Popen, DEVNULL
from distutils.spawn import find_executable
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO, SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import UnknownSeq, Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter


def create_phylogeny_wrapper(wd, **kwargs):

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
# mafft
if find_executable('mafft') is not None:
    pass
else: 
    print('You do not have MAFFT installed in your path. Please fix this.')
    sys.exit()
# fasttree, later check if phytool is raxml and check for correct raxml argument
if find_executable('FastTree') is not None:
    pass
else: 
    print('You do not have FastTree installed in your path. Please fix this.')
    sys.exit()
