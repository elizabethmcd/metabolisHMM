from setuptools import setup, find_packages
import os

with open('requirements.txt') as f:
    requirements= f.readlines()

setup(name='metabolisHMM',
version=2.22,
description='Constructing phylogenies and performing functional annotations with HMM markers',
url='https://github.com/elizabethmcd/metabolisHMM',
author='Elizabeth McDaniel',
author_email='emcdaniel@wisc.edu',
license='GNU',
packages=find_packages(),
scripts=['bin/create-genome-phylogeny','bin/search-custom-markers','bin/single-marker-phylogeny', 'bin/summarize-metabolism'],
install_requires=requirements,
zip_safe=False)
