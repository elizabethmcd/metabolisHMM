# metabolisHMM

A tool for summarizing metabolic capabilities and creating phylogenies of genomes based on manually curated markers.

## Dependencies and Setup 

- python (Biopython)
- HMMer
- MUSCLE/Mafft
- [catfasta2phyml](https://github.com/nylander/catfasta2phyml)
- FastTree/RaxML

To run the pipeline, clone the repository with `git clone https://github.com/elizabethmcd/metabolisHMM`. This will provide you with the necessary scripts and HMM markers to run the analysis on a given set of genomes. If your filenames have dashes `-` in them, rename them for example with underscores instead, with `rename 's/-/_/g' *`. You will need all of the above programs in your path, including the `catfasta2phyml` perl script.

## Metabolic Summaries

To get metabolic summaries of your genomes, run `python summarize-metabolism.py`. 
```
usage: python summarize-metabolism.py --genome_dir --output
  --genome_dir Directory where genomes to be screened are held
  --output Name of summary results file
```

This will create a summary table of each metabolic marker and how many hits were found above the given threshold in all of your genomes. 

If you want to run the summary on a specific set of markers outside of those included here (such as a specific pathway - e.g. Wood Ljungdhal as provided), use `python search-custom-markers.py`. 
```
usage: python search-custom-markers.py --genome_dir --markers_dir --output
  --genome_dir Directory where genomes to be screened are held
  --markers_dir Directory where custom markers are held
  --output Name of results file
```

This will place the columns in alphabetical order of markers, and not necessarily the order of the pathway.

## Genome Phylogenies

To create a genome phylogeny on the given genomes based on ribosomal protein markers, run `python create-genome-phylogeny.py`. 

```
usage: python create-genome-phylogeny.py --genome_dir --domain --phylogeny --threads
Creates full archaeal/bacterial genome phylogenies based off specific ribosomal protein markers
arguments:
  --genome_dir Directory where genomes to be screened are held
  --domain archaea, bacteria
  --phylogeny fastree, raxml
  --threads #threads for performing alignments and calculating phylogeny
```

## Phylogeny of a single marker

Say you are interested in the phylogenetic distribution of a particular marker (nifA for example) amongst your genomes. Use the script `single-marker-phylogeny.py` and the path to your marker of interest. You can put it in the `metabolic-markers` directory to be run with others for metabolic summaries, or have it in a different location. The usage is: 

```
usage: python single-marker-phylogeny.py --genome_dir --marker --list --phylogeny --threads
Creates phylogeny of a single marker against given set of genomes
arguments:
  --genome_dir Directory where genomes to be screened are held
  --marker directory/name of specific marker to be screened
  --list name of list with marker hits
  --phylogeny fastree, raxml
  --threads #threads for calculating phylogeny
```

## Caveats 

This tool is still under active development, and may contain bugs here and there, so you will want to perform some sanity-checks manually. All HMMer results are in the `out` folder and split by genome. For metabolic marker genes, hits are recordered if they are above the "TC" threshold that we manually defined. For ribosomal marker genes for making genome phylogenies, the best hit is taken so that for every genome and ribosomal marker, only one hit is taken into account. For MAGs with multiple hits, you may want to check if a certain hit falls on the end of a contig and is not truly the "best" hit.  

## References 

Some of the markers were made publicly available from Karthik Anantharaman's study on [Thousands of microbial genomes shed light on interconnected biochemical processes in an aquifer system](https://www.nature.com/articles/ncomms13219) found [here](https://github.com/kanantharaman/metabolic-hmms). Certainly using HMM profiles to characterize metabolic capabilities predates this publication, but if you wanted you could cite this paper as a proof of practice. 

1. Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system. K. Anantharaman et al. _Nature Communications_. 2016. 

Additionally if you add a new marker to the database, and your analysis heavily focusses on that marker (e.g. nitrogen cycling), you may want to cite some reference about that marker protein. 

The ribosomal protein markers were also used for the expanded view of the tree of life from Hug et al.:

2. A new view of the tree of life. L. Hug et al. _Nature Microbiology_. 2016.

All credit for the concatenation script goes to John Nylander, as it works very nicely and quickly for these purposes. 

The code and released HMM profiles are free to use for any purposes. Feel free to submit an issue with any questions or something I glaringly did wrong! 