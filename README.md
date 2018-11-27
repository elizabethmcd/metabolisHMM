# metabolisHMM

A tool for summarizing the metabolic capabilities and creating phylogenies of a set of genomes based on manually curated marker genes. 

## Dependencies and Setup 

- python
  - BioPython
- HMMer
- Muscle

To run the pipeline, clone the repository with `git clone https://github.com/elizabethmcd/metabolisHMM`. This will provide you with the necessary scripts and HMM markers to run the analysis on a given set of genomes. Place the protein files ending in the extension `.faa` in a folder named `genomes`. 

To get metabolic summaries of your genomes, run `python run-metabolic-markers.py`. This will create a summary table of each metabolic marker and how many hits were found above the given threshold in all of your genomes. To create a genome phylogeny on the given genomes based on ribosomal protein markers, run `python run-ribosomal-markers.py`. 

## Caveats 

## References 