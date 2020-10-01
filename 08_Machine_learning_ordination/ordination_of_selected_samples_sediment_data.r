#!/usr/bin/env Rscript

# unsupervised ordination of sediment parameters
library(speedyseq)
library(ggplot2)
library(ranger)
library(vegan)

# set paths
projectPath <- "/data/projects/2019/tnt"
resultPath <- file.path(projectPath, "ML/results")
phyloPath <- file.path(projectPath,"phyloseq_output")
# read in phyloseq object
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))
# set parameters
num.trees <- 10000
mtry_factor <- 1
ML <- "ranger"
nodename <- Sys.info()[["nodename"]]

# get sample data table and select columns for ordination
correlation_cols <- c("Hg_microg_kg", "Pb206_207_ratio", "P_percent", "Ca_percent",
  "Fe_percent", "Sc45_ppm","V51_ppm", "Cr52_ppm","Mn55_ppm","Co59_ppm","Ni60_ppm","Cu63_ppm","Zn66_ppm","As75_ppm", 
  "Sr86_ppm","Zr90_ppm","Mo95_ppm","Ag107_ppm", "Cd111_ppm", "Sn118_ppm", "Sb121_ppm", "Cs133_ppm", "Ba137_ppm", "W186_ppm", 
  "Tl205_ppm", "Pb207_ppm", "Bi209_ppm", "Th232_ppm", "U238_ppm", 
  "TIC", "TN", "TC", "TS", "TOC", 
  "Microm_001_mean", "Microm_63_mean","Microm_125_mean", "Microm_250_mean", "Microm_500_mean", "Microm_1000_mean",
  "TNT", "ADNT_2", "ADNT_4", "DANT_2.4", "DANT_2.6", "DNT_2.4", "DNT_2.6", "DNB", "TNB", "HMX", "RDX", "Tetryl")
sample_data <- as.data.frame(unclass(sample_data(ps_ML)))
row.names(sample_data) <- row.names(sample_data(ps_ML))
sample_data2 <- sample_data[, names(sample_data) %in% correlation_cols]

# prepare data for unsupervised random forest and get mtry number
unsupervised_table <- sample_data2
all_vars <- ncol(unsupervised_table)
for_mtry <- ifelse((sqrt(all_vars) * mtry_factor) < all_vars,
  sqrt(all_vars) * mtry_factor, all_vars)
# generate noised data set with shuffled columns
synth_data <- as.data.frame(lapply(as.data.frame(unsupervised_table), function(x) {
  sample(x, length(x), replace = TRUE)
}))
# combine original and noised data set
combined_data <- rbind(data.frame(y = 0, unsupervised_table), 
             data.frame(y = 1, synth_data))
combined_data$y <- factor(combined_data$y)

# Run unsupervised Random Forest and extract proximity matriy (helper function)
forest <- ranger(y ~ ., combined_data, keep.inbag = TRUE, num.trees = num.trees, mtry = for_mtry)
prox <- extract_proximity_oob(forest, combined_data)[1:nrow(unsupervised_table), 1:nrow(unsupervised_table)]
row.names(prox) <- row.names(unsupervised_table)

# Generate dist matrix and % of variation for x and y axis
distance_matrix <- dist(1 - prox)
PCA_object <- cmdscale(distance_matrix, eig = TRUE, x.ret = TRUE)

PCA_variation <- round(PCA_object$eig/sum(PCA_object$eig) * 100, 1)
PCA_values <- scores(PCA_object)

# Combine to plot ready data frame
PCA_data <- data.frame(
  Sample = rownames(PCA_values),
  X = PCA_values[, 1],
  Y = PCA_values[, 2])

# merge with sample data
RF_PCA_meta <- merge(PCA_data, sample_data, 
  by.x = "Sample", by.y = "row.names", all.x = TRUE)
RF_PCA_meta$Target <- cut(RF_PCA_meta$TNT, breaks = c(-Inf, 0, Inf), labels = c("absent", "present"))
RF_PCA_meta$Nodename <- nodename

write.table(RF_PCA_meta, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", nodename, "_URF_PCA_sediment_full.tsv")))
write.table(PCA_variation, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", nodename, "_URF_PCA_sediment_full_variation.tsv")))
