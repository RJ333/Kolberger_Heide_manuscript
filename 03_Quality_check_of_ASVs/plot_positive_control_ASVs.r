#!/usr/bin/env Rscript

# This script identifies the positive controls used in PCR for all plates
library(speedyseq)
library(ggplot2)
library(rlang)

# select project directory,
# the phyloseq object and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
phyloPath <- file.path(projectPath,"phyloseq_output")
plotPath <- file.path(projectPath, "plots")
# read in phyloseq object
ps_full <- readRDS(file.path(phyloPath, "ps_full.RDS"))

# identify which strains were used as which positive controls
positives <- subset_samples(ps_full, Library_purpose == "positive_pcr_control")
positive_controls <- filter_taxa(positives, function (x) {sum(x > 50) >= 1}, prune = TRUE)      
positive_controls_melt <- psmelt(positive_controls)

# plotting requires helper function "position_jitter_stack"
ggplot(positive_controls_melt, aes(x = Sample, y = Abundance, fill = Genus, colour = Genus)) +
  geom_bar(stat = "identity", colour = "black") + 
  geom_text(aes(label = OTU), position = position_jitter_stack(
      vjust =0.5, jitter.height = 1000, jitter.width =  0.4, offset = 1), 
    size = 2, alpha = 0.7, colour = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
