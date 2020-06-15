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
    subparsers = parser.add_subparsers(help='Desired workflow',dest='operation')

    parent_parser = argparse.ArgumentParser(add_help=False)
    flags = parent_parser.add_argument_group('SYSTEM PARAMETERS')
    flags.add_argument('-p','--processors',help='threads',default=2,type=int)
    flags.add_argument('-d','--debug',help='make extra debugging output',default=False,
                        action= "store_true")
    flags.add_argument("-h", "--help", action="help", help="show this help message and exit")

 # Parser for each workflow
    
    # Genome phylogeny parser 
    genome_phylogeny_parent = argparse.ArgumentParser(add_help=False)
    required_phylogeny_flags = genome_phylogeny_parent.add_argument_group('REQUIRED ARGUMENTS')
    required_phylogeny_flags.add_argument("-i","--input", help="Directory where genomes to be screend are held")
    required_phylogeny_flags.add_argument("-o","--output", help="Directory to store results and intermediate files")
    required_phylogeny_flags.add_argument("-dom","--domain", help="Options: Archaea or Bacteria (each requires separate curated markers and thus must be constructed separately")
    required_phylogeny_flags.add_argument("-phy","--phylogeny", help="Options: fastree or raxml")
    metadata_phylogeny_flags = genome_phylogeny_parent.add_argument_group('METADATA ARGUMENTS')
    metadata_phylogeny_flags.add_argument("-loc","--loci", help="Alert the user if a certain genome contains less than X number of ribosomal protein loci. Default=12", default=12)
    metadata_phylogeny_flags.add_argument("-n", "--names", help="CSV formatted file of filenames corresponding to taxonomical or group names for each genome")
    metadata_phylogeny_flags.add_argument("-itol", "--itol_file", help="Output ITOL formatted metadata file. Default = itol_metadata.txt", default="itol_metadata.txt")

# Arguments for each operation

    # Genome phylogeny operation
    create_phylogeny_parser = subparsers.add_parser("create_phylogeny",formatter_class=SmartFormatter,parents=[parent_parser,genome_phylogeny_parent], add_help=False)

    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        printHelp()
        sys.exit(1)
    else:
        return parser.parse_args(args)
