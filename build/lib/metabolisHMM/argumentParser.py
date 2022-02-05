#!/usr/bin/env python3

###############################################
# metabolisHMM argument parse for workflows
###############################################

__author__="Elizabeth McDaniel"
__license__="GPL"
__email__="elizabethmcd93@gmail.com"


import argparse
import os
import sys

import metabolisHMM
from metabolisHMM.controller import Controller

def version():
    versionFile = open(os.path.join(metabolisHMM.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()

class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

# Top level help menu

def printHelp():
    print('')
    print('         metabolisHMM v' + VERSION + '    ')
    print('''\

  Elizabeth McDaniel. GPL 3.0 License. McMahon Lab, UW-Madison. 2020
  
  Choose one of the workflows below for more detailed help. See https://github.com/elizabethmcd/metabolisHMM/wiki for documentation
  
  Example: metabolisHMM custom_search -h

  Workflows:
        create_phylogeny            -> Create genome phylogeny of curated ribosomal protein markers
        single_marker_phylogeny     -> Create a phylogeny of a single marker
        summarize_metabolism        -> Summarize metabolic capabilities with curated markers
        custom_search               -> Search for custom sets of markers among genomes

  ** Genome names can only contain underscores "_" and not hypens "-" **

    ''')

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class = SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired workflow')


    create_phylogeny_parser = subparsers.add_parser("create_phylogeny",formatter_class=SmartFormatter, add_help=False)
    create_phylogeny_parser.add_argument("-i","--input",help="Directory where genomes to be screened are held")

 # Parser for each workflow
    
    # Genome phylogeny parser 
    genome_phylogeny_parent = argparse.ArgumentParser(add_help=False)
    phyflags = genome_phylogeny_parent.add_argument_group('REQUIRED ARGUMENTS')
    phyflags.add_argument("-i","--input", help="Directory where genomes to be screend are held")
    phyflags.add_argument("-o","--output", help="Directory to store results and intermediate files")
    phyflags.add_argument("-dom","--domain", help="Options: Archaea or Bacteria (each requires separate curated markers and thus must be constructed separately")
    phyflags.add_argument("-phy","--phylogeny", help="Options: fastree or raxml")

# Arguments for each operation

    # Genome phylogeny operation
    
    # Genome phylogeny metadata flags
    metadata = create_phylogeny_parser.add_argument_group("METADATA ARGUMENTS")
    metadata.add_argument("-loc","--loci", help="Alert the user if a certain genome contains less than X number of ribosomal protein loci. Default=12", default=12)
    metadata.add_argument("-n", "--names", help="CSV formatted file of filenames corresponding to taxonomical or group names for each genome")
    metadata.add_argument("-itol", "--itol_file", help="Output ITOL formatted metadata file. Default = itol_metadata.txt", default="itol_metadata.txt")
    
    
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        printHelp()
        sys.exit(1)
    else:
        return parser.parse_args()
