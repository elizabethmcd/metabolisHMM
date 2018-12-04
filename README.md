# metabolisHMM

A tool for summarizing the metabolic capabilities and creating phylogenies of a set of genomes based on manually curated marker genes. 

**Note: At the present time, please email me at emcdaniel@wisc.edu for access to metabolic markers, as this is heavly under development and not ready for relase.**

## Dependencies and Setup 

- python
  - BioPython
  - pandas
- HMMer
- Muscle

To run the pipeline, clone the repository with `git clone https://github.com/elizabethmcd/metabolisHMM`. This will provide you with the necessary scripts and HMM markers to run the analysis on a given set of genomes. Place the protein files ending in the extension `.faa` in a folder named `genomes`. 

To get metabolic summaries of your genomes, run `python run-metabolic-markers.py`. This will create a summary table of each metabolic marker and how many hits were found above the given threshold in all of your genomes. To create a genome phylogeny on the given genomes based on ribosomal protein markers, run `python run-ribosomal-markers.py`. Metabolic characterization summaries are provided in the `results` folder, and columns are ordered by functions *which will be added soon*. 

## Caveats 

This tool is still under active development, and may contain bugs here and there, so you will want to perform some sanity-checks manually. All HMMer results are in the `out` folder and split by genome. For metabolic marker genes, hits are recordered if they are above the "TC" threshold that we manually defined. **Known bug: 4 of the HMMs don't have TC thresholds yet, and therefore do not run with the others, yet.**. For ribosomal marker genes for making genome phylogenies, the best hit is taken so that for every genome and ribosomal marker, only one hit is taken into account. For MAGs with multiple hits, you may want to check if a certain hit falls on the end of a contig and is not truly the "best" hit. **As of 2018-12-03, the full `run-ribosomal-markers.py` pipeline is not finished beyond running HMMs.**. 

## References 

Some of the markers were made publicly available from Karthik Anantharaman's study on [Thousands of microbial genomes shed light on interconnected biochemical processes in an aquifer system](https://www.nature.com/articles/ncomms13219) found [here](https://github.com/kanantharaman/metabolic-hmms). Certainly using HMM profiles to characterize metabolic capabilities predates this publication, but if you wanted you could cite this paper as a proof of practice. 

1. Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system. K. Anantharaman et al. _Nature Communications_, 2016. 

Additionally if you add a new marker to the database, and your analysis heavily focusses on that marker (e.g. nitrogen cycling), you may want to cite some reference about that marker protein. 

The code and HMM profiles are free to use for any purposes. Feel free to submit an issue with any questions or something I glaringly did wrong! 