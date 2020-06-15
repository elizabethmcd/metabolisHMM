#!/usr/bin/env python3

###############################################
# metabolisHMM program entry point to different workflows
###############################################

__author__="Elizabeth McDaniel"
__license__="GPL"
__email__="elizabethmcd93@gmail.com"

import argparse
import logging
import os
import sys

import metabolisHMM
import metabolisHMM.create_phylogeny

def version():
    versionFile = open(os.path.join(metabolisHMM.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()

class Controller():
    def __init__(self):
        self.logger = logging.getLogger()

    def create_phylogeny_operation(self, **kwargs):
        logging.debug("Creating genome phylogeny...")
        metabolisHMM.create_phylogeny.create_phylogeny_wrapper(**kwargs)
        logging.debug("Finished creating genome phylogeny!")
    
    def custom_search_operation(self, **kwargs):
        logging.debug("Starting search of custom markers...")
        #
        logging.debug("Finished custom search!")
    
    def single_marker_phylogeny_operation(self, **kwargs):
        logging.debug("Starting phylogeny of single marker...")
        #
        logging.debug("Finished single marker phylogeny!")
    
    def summarize_metabolism_operation(self, **kwargs):
        logging.debug("Summarizing metabolism of curated markers...")
        #
        logging.debug("Finished summarizing metbolic markers!")
    
    def setup_logger(self,loc):

        # set up logging everything to file
        logging.basicConfig(level=logging.DEBUG,
                           format='%(asctime)s %(levelname)-8s %(message)s',
                           datefmt='%m-%d %H:%M',
                           filename=loc)

        # set up logging of INFO or higher to sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        logging.debug("***Logger started at {0}***".format(loc))
        logging.debug("Command to run metabolisHMM was: {0}\n".format(' '.join(sys.argv)))
        logging.debug("metabolisHMM version {0} was run \n".format(VERSION))   
    
    def parseArguments(self, args):
        # Correct workflow

        if args.operation == "create_phylogeny":
            self.create_phylogeny_operation(**vars(args))
        if args.operation == "custom_search":
            self.custom_search_operation(**vars(args))
        if args.operation == "single_marker_phylogeny":
            self.single_marker_phylogeny_operation(**vars(args))
        if args.operation == "summarize_metabolism":
            self.summarize_metabolism_operation(**vars(args))
    
    def loadDefaultArgs(self):
        pass