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
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

# Arguments passed from summarize-metabolism.py with metabolic stats output, and input from user for metadata and aggregating option to visualize

userprefs <- commandArgs(trailingOnly = TRUE)
markers.file.path <- userprefs[1]
metadata.file.path <- userprefs[2]
aggregate.option <- userprefs[3]
heatmap.file.path <- userprefs[4]
heatmap.option <- userprefs[5]

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

make.heatmap <- function(stats_table, file.path, heatmap.option, directory){
  option <- heatmap.option
  if(option == 'CUSTOM'){
    heatmap.file.path <- file.path
    figureOut <- file.path(directory,heatmap.file.path)
    table_melted <- melt(stats_table, id.vars="group")
    table_melted <- table_melted %>% mutate(group = factor(group),group = factor(group, levels = rev(levels(group))))
    plot <- table_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black') + scale_fill_viridis_c(alpha=1,begin=0,end=1,direction=-1) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    plot_formatted <- plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=85, hjust=0), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + labs(x=NULL,y=NULL) + guides(fill = guide_colorbar(nbin = 10)) + scale_y_discrete(expand=c(0,0))
    ggsave(file=figureOut, plot_formatted, height=15, width=55, units=c("cm"))
  }
  if(heatmap.option == 'CURATED'){
    heatmap.file.path <- file.path
    figureOut <- file.path(directory,heatmap.file.path)
    carbon = cbind(stats_table$group, presence_absence[,c(1:26)])
    colnames(carbon) <- c("group", "MtOH dehy", "madA", "madB", "fdh", "sfh", "sgdh", "smdh", "fae", "fmtF", "mtmc", "fdhA", "fdhB", "fdhC", "coxL", "coxM", "coxS", "rubisco I", "rubisco II", "rubisco III", "rubisco II/III", "rubisco IV", "codhC", "codhD", "codh catalytic", "aclA", "aclB")
    nitrogen = cbind(stats_table$group, presence_absence[,c(27:47)])
    colnames(nitrogen) <- c("group", "nifA", "nifB", "nifH", "nxrA", "nxrB", "napA", "napB", "narG", "narH", "nrfA", "nrfH", "nirB", "nirD", "nirK", "nirS", "norB", "norC", "nosD", "nosZ", "hzoA", "hzsA")
    sulfur = cbind(stats_table$group, presence_absence[,c(48:59)])
    colnames(sulfur) <- c("group", "fccB", "sqr", "sdo", "aprA", "sat", "dsrA", "dsrB", "dsrD", "phsA", "soxB", "soxC", "soxY")
    oxygen = cbind(stats_table$group, presence_absence[,c(60:70)])
    colnames(oxygen) <- c("group", "coxA", "coxB", "ccoN", "ccoO", "ccoP", "cyoA", "cyoD", "cyoE", "cydA", "cydB", "qoxA")
    hydrogen = cbind(stats_table$group, presence_absence[,c(71:80)])
    colnames(hydrogen) <- c("group", "FeFe I", "FeFe II", "Group I", "Group IIA", "Group IIB", "Group IIIA", "Group IIIB", "Group IIIC", "Group IIID", "Group IV")
    # individual melted dfs and plots
    # carbon
    carbon_melted <- melt(carbon, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group))))
    carbon_plot <- carbon_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() +  scale_fill_gradient(low="gray92", high="brown4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    carbon_plot_formatted <- carbon_plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + labs(x="Carbon", y=NULL) +  scale_y_discrete(expand=c(0,0))
    carbon_clean <- carbon_plot_formatted + theme(legend.position="none")
    # nitrogen
    nitrogen_melted <- melt(nitrogen, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group)))) 
    nitrogen_plot <- nitrogen_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black', size=0.5,aes(width=1, height=1)) + coord_fixed() + scale_fill_gradient(low="gray92", high="dodgerblue4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    nitrogen_plot_formatted <- nitrogen_plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Nitrogen", y=NULL) +  scale_y_discrete(expand=c(0,0))
    nitrogen_clean <- nitrogen_plot_formatted + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    # sulfur
    sulfur_melted <- melt(sulfur, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group)))) 
    sulfur_plot <- sulfur_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() + scale_fill_gradient(low="gray92", high="darkorchid4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    sulfur_plot_formatted <- sulfur_plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Sulfur", y=NULL) +  scale_y_discrete(expand=c(0,0))
    sulfur_clean <- sulfur_plot_formatted + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    # oxygen
    oxygen_melted <- melt(oxygen, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group)))) 
    oxygen_plot <- oxygen_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() + scale_fill_gradient(low="gray92", high="turquoise4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    oxygen_plot_formatted <- oxygen_plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Oxygen", y=NULL) +  scale_y_discrete(expand=c(0,0))
    oxygen_clean <- oxygen_plot_formatted + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
    # hydrogen
    hydrogen_melted <- melt(hydrogen, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group)))) 
    hydrogen_plot <- hydrogen_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black', size=0.5,aes(width=1, height=1)) + coord_fixed() + scale_fill_gradient(low="gray92", high="darkorange3") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.grid = element_blank(), panel.border = element_blank())
    hydrogen_plot_formatted <- hydrogen_plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Hydrogen", y=NULL) +  scale_y_discrete(expand=c(0,0))
    hydrogen_clean <- hydrogen_plot_formatted + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank())
    # stitch together
    carbon_tmp <- ggplot_build(carbon_clean)
    nitrogen_tmp <- ggplot_build(nitrogen_clean)
    sulfur_tmp <- ggplot_build(sulfur_clean)
    oxygen_tmp <- ggplot_build(oxygen_clean)
    hydrogen_tmp <- ggplot_build(hydrogen_clean)
    g1 <- ggplot_gtable(carbon_tmp) ; g2 <- ggplot_gtable(nitrogen_tmp) ; g3 <- ggplot_gtable(sulfur_tmp) ; g4 <- ggplot_gtable(oxygen_tmp) ; g5 <- ggplot_gtable(hydrogen_tmp)
    n1 <- length(carbon_tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
    n2 <- length(nitrogen_tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
    n3 <- length(sulfur_tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
    n4 <- length(oxygen_tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
    n5 <- length(hydrogen_tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
    g <- cbind(g1, g2, g3, g4, g5, size="first")
    ggsave(file=figureOut, plot=g, width=40, height=40, units=c("cm"))
  }
}

# call functions
markers <- import.markers(File = markers.file.path)
metadata <- import.metadata(File = metadata.file.path)
out <- get.output.directory(file.path=markers.file.path)
stats_table <- merge.markers.metadata(markers, metadata, aggregate.option, outdir=out)
heatmap <- make.heatmap(stats_table, file.path=heatmap.file.path, heatmap.option, directory=out )
