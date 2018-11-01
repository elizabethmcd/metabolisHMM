#! /usr/bin/env/ python 

import os 

def version():
    """Read program version from file."""
    version_file = open(os.path.join_path_[0], "VERSION")
    return version_file.readline.strip()