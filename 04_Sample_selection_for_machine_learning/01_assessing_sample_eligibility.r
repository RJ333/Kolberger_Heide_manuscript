#!/usr/bin/env Rscript

# Select and prepare the libraries which can be used for Machine Learning
# without excessive risk of confounding variables
library(ggplot2)
library(dplyr)
library(speedyseq)

# select project directory, phyloseq path and the save path
projectPath <- "/data/projects/2019/tnt"
savePath <- file.path(projectPath, "ML/input_tables")
phyloPath <- file.path(projectPath,"phyloseq_output")
# read in phyloseq object and extract sample data
ps_final <- readRDS(file.path(phyloPath, "ps_final.RDS"))
samples <- as.data.frame(unclass(sample_data(ps_final)))
samples$TNT_cut <- cut(samples$TNT, breaks = c(-Inf, 0, Inf), labels = c("absent", "present"))

# identify samples too similar based on criteria in manuscript methods section
samples %>% 
  group_by(Primerset, Udemm, Nucleic_acid, Station_ID, Collection, Cruise_ID, Biological_replicate, Experiment, Sample_Type, Area) %>% 
  select(TNT_cut) %>%
  distinct() %>%
  ungroup() %>%
  group_by(Sample_Type, Cruise_ID, Experiment, Area, Primerset, Nucleic_acid, Collection, Station_ID) %>% # Station ID to see which samples
  count(TNT_cut)
