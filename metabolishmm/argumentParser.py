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

import metabolishmm
from metabolishmm.controller import Controller

def version():
    versionFile = open(os.path.join(metabolishmm.__path__[0], 'VERSION'))
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
    flags = parent_parser.add_argument_group('System paramters')
    flags.add_argument('-p','--processors',help='threads',default=2,type=int)
    flags.add_argument('-d','--debug',help='make extra debugging output',default=False,
                        action= "store_true")
    flags.add_argument("-h", "--help", action="help", help="show this help message and exit")

 # Parser for each workflow
    
    genome_parent = argparse.ArgumentParser(add_help=False)
    GenomeFlags = genome_parent.add_argument_group('Create a genome phylogeny using curated ribosomal protein markers')


