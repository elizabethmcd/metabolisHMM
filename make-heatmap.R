#! /usr/bin/env R

#####################################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# make-heatmap : visualize output of metabolic summaries (given or custom set) in a heatmap
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# November 2018
# This program is free software under the GNU General Public License version 3.0
#####################################################################
suppressMessages(library(tidyverse))
message = FALSE

# Arguments passed from summarize-metabolism.py with metabolic stats output, and input from user for metadata and aggregating option to visualize

userprefs <- commandArgs(trailingOnly = TRUE)
markers.file.path <- userprefs[1]
metadata.file.path <- userprefs[2]
aggregate.option <- userprefs[3]

import.markers <- function(File){
  markers.file.path <- File
  markers <- read.csv(file <- markers.file.path, header = TRUE)
  colnames(markers)[1] <- "genome"
  return(markers)
}

import.metadata <- function(File){
  metadata.file.path <- File
  metadata <- read.csv(file <- metadata.file.path, header = FALSE)
  colnames(metadata) <- c("genome", "group")
  return(metadata)
}

merge.markers.metadata <- function(markers, metadata, aggregate.option){
  merged_table <- left_join(metadata, markers, by="genome")
  presence_absence = as.data.frame(lapply(merged_table[3:ncol(merged_table)], function(x) ifelse(x>1, 1, x)))
  presence_absence$group <- metadata$group
  exclude <- ncol(presence_absence) - 1
  sums <- as.data.frame(colMeans)
  if(aggregate.option == "OFF"){
    # add column for marker avg % in all genomes
    return(presence_absence)
  } else {
    exclude <- ncol(presence_absence) - 1
    agg <- aggregate(presence_absence[1:exclude], list(presence_absence$group), mean)
    colnames(agg)[1] <- "group"
    return(agg)
  }
}

# call functions
markers <- import.markers(File = markers.file.path)
metadata <- import.metadata(File = metadata.file.path)
result <- merge.markers.metadata(markers, metadata, aggregate.option)
print(result)