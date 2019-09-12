#! /usr/bin/env R

#####################################################################
# metabolisHMM - A tool for exploring and visualizing the distribution and evolutionary histories of metabolic markers
# make-heatmap : visualize output of metabolic summaries (given or custom set) in a heatmap
# Written by Elizabeth McDaniel emcdaniel@wisc.edu
# September 2019
# This program is free software under the GNU General Public License version 3.0
#####################################################################
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))

# Arguments passed from summarize-metabolism.py with metabolic stats output, and input from user for metadata and aggregating option to visualize

userprefs <- commandArgs(trailingOnly = TRUE)
markers.file.path <- userprefs[1]
metadata.file.path <- userprefs[2]
aggregate.option <- userprefs[3]
heatmap.file.path <- userprefs[4]

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

get.output.directory <- function(file.path){
  markers.file.path <- file.path
  dir <- dirname(markers.file.path)
  return(dir)
}


merge.markers.metadata <- function(markers, metadata, aggregate.option, outdir){
  merged_table <- left_join(metadata, markers, by="genome")
  presence_absence = as.data.frame(lapply(merged_table[3:ncol(merged_table)], function(x) ifelse(x>1, 1, x)))
  presence_absence$group <- metadata$group
  if(aggregate.option == "OFF"){
    out <- outdir
    outname <- file.path(out,"cleaned-metabolic-table.csv")
    write.csv(file=outname, presence_absence, row.names=FALSE, quote=FALSE)
    return(presence_absence)
  } else {
    exclude <- ncol(presence_absence) - 1
    agg <- aggregate(presence_absence[1:exclude], list(presence_absence$group), mean)
    colnames(agg)[1] <- "group"
    out <- outdir
    outname <- file.path(out,"cleaned-aggregated-metabolic-table.csv")
    write.csv(file=outname,agg, row.names=FALSE, quote=FALSE)
    return(agg)
  }
}

make.heatmap <- function(stats_table, file.path, directory){
  heatmap.file.path <- file.path
  figureOut <- file.path(directory,heatmap.file.path)
  table_melted <- melt(stats_table, id.vars="group")
  plot <- table_melted %>% ggplot(aes(x=variable, y=group, fill=value)) + geom_tile(color='black') + scale_fill_viridis_c(alpha=1,begin=0,end=1,direction=-1) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())
  plot_formatted <- plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10)) + scale_y_discrete(expand=c(0,0))
  ggsave(file=figureOut, plot_formatted, height=15, width=55, units=c("cm"))
}

# call functions
markers <- import.markers(File = markers.file.path)
metadata <- import.metadata(File = metadata.file.path)
out <- get.output.directory(file.path=markers.file.path)
result <- merge.markers.metadata(markers, metadata, aggregate.option, outdir=out)
make.heatmap(stats_table = result, file.path = heatmap.file.path, directory = out)
