<img align="right" width="400" src="https://github.com/elizabethmcd/metabolisHMM/blob/master/data/metabolisHMM-logo.png">

<br>
<br>
<br>

---

# metabolisHMM: a tool for summarizing metabolic distributions and evolutionary relationships

[metabolisHMM](https://github.com/elizabethmcd/metabolisHMM/wiki) is a pipeline for visualizing the distribution and evolutionary relationships of specific metabolic markers using Hidden Markov Model (HMM) profiles. Annotation through HMMs is becoming increasingly popular over BLAST-based methods, however methods for rapidly visualizing and summarizing these results is lacking. This software automates the process of searching any set of metabolic markers that have an HMM built against a set of genomes, and outputs phylogenies, summary statistics, and heatmap distributions. The metabolisHMM software is primarily written in python for performing pipeline steps and parsing results, with R visualization steps added for producing heatmaps. 

See the [wiki](https://github.com/elizabethmcd/metabolisHMM/wiki) for installation and usage instructions. 

## Dependencies and Setup 

- [python (>3.6) (Biopython, pandas)](https://www.anaconda.com/)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [HMMER](http://hmmer.org/)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [RaxML](https://cme.h-its.org/exelixis/software.html)
- [TrimAl](https://github.com/scapella/trimal)
- [R (tidyverse, reshape2)](https://cran.r-project.org/)

To run the pipeline, clone the repository with `git clone https://github.com/elizabethmcd/metabolisHMM`. This will provide you with the necessary scripts and HMM markers to run the analysis on a given set of genomes. If your filenames have dashes `-` in them, rename them for example with underscores instead, with `rename 's/-/_/g' *`.

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
  --input Directory where genomes to be screened are held
  --output Directory to hold intermediate files and results
  --domain archaea, bacteria, all
  --phylogeny fastree, raxml
  --threads #threads for performing alignments and calculating phylogeny with raxml, doesn't matter for fastree
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

## References 

Some of the markers were made publicly available from Karthik Anantharaman's study on [Thousands of microbial genomes shed light on interconnected biochemical processes in an aquifer system](https://www.nature.com/articles/ncomms13219) found [here](https://github.com/kanantharaman/metabolic-hmms). Certainly using HMM profiles to characterize metabolic capabilities predates this publication, but if you wanted you could cite this paper as a proof of practice. 

1. Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system. K. Anantharaman et al. _Nature Communications_. 2016. 

Additionally if you add a new marker to the database, and your analysis heavily focusses on that marker (e.g. nitrogen cycling), you may want to cite some reference about that marker protein. 

The ribosomal protein markers were also used for the expanded view of the tree of life from Hug et al.:

2. A new view of the tree of life. L. Hug et al. _Nature Microbiology_. 2016.

The code and released HMM profiles are free to use for any purposes. Feel free to submit an issue with any questions or something I glaringly did wrong! 