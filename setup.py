from setuptools import setup, find_packages
import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'metabolisHMM', 'VERSION'))
    return versionFile.read().strip()

setup(name='metabolisHMM',
version=version(),
description='Constructing phylogenies and performing functional annotations with curated and custom HMM markers',
url='https://github.com/elizabethmcd/metabolisHMM',
author='Elizabeth McDaniel',
author_email='elizabethmcd93@gmail.com',
license='GNU',
package_data={'metabolisHMM': ['VERSION']},
packages=['metabolisHMM'],
scripts=['bin/metabolisHMM'],
install_requires=[
    'numpy',
    'pandas',
    'biopython',
    'seaborn',
    'matplotlib'
],
zip_safe=False)
