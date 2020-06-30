#!/usr/bin/env python3

import os
import sys
import logging


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('timestamp')
        self.timeKeeper = TimeKeeper()

    def create_phylogeny(self, options):
    
    def single_marker(self, options):
    
    def custom_search(self, options):

    def summarize_metabolism(self, options):
    
    
    def parseOptions(self, options):

        if(options.subparser_name == 'create-phylogeny'):
            self.create_phylogeny(options)
        elif(options.subparser_name == 'single-marker'):
            self.single_marker(options)
        elif(options.subparser_name == 'custom-search'):
            self.custom_search(options)
        elif(options.subparser_name == 'summarize-metabolism'):
            self.summarize_metabolism(options)
        else:
            self.logger.error("Unknown metabolisHMM command: " + options.subparser_name + '\n')
            sys.exit(1)
        
        return 0