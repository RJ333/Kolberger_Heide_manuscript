#!/usr/bin/env Rscript
# This script prepares the data to predict the presence of TNT only using sediment data without 16S communities
library(speedyseq)
library(phyloseq2ML)
library(futile.logger)
flog.threshold(INFO)

projectPath <- "/data/projects/2019/tnt"
savePath <- file.path(projectPath, "ML/input_tables")
phyloPath <- file.path(projectPath,"phyloseq_output")
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))

# settings for data preparation
ML <- "ranger"
mode <- "classification"
analysis_ID <- "ave_run10"
split_ratios <- 0.75
my_seed <- readr::parse_number(Sys.info()[["nodename"]])  # seed for chandler-1 is 1
copies <- 0
noise_factor  <- 0
log_time_mem <- TRUE
desired_response_vars <- "TNT"

# if desired, subset the correlation cols by these values
best_correlation_cols <- c("As75_ppm", "Microm_63_mean", "Microm_250_mean", 
  "Microm_500_mean", "V51_ppm", "Zn66_ppm", "TN", "Co59_ppm", "Fe_percent")

# get sample data table and select columns
correlation_cols <- c("Hg_microg_kg", "Pb206_207_ratio", "P_percent", "Ca_percent",
  "Fe_percent", "Sc45_ppm","V51_ppm", "Cr52_ppm","Mn55_ppm","Co59_ppm","Ni60_ppm","Cu63_ppm","Zn66_ppm","As75_ppm", 
  "Sr86_ppm","Zr90_ppm","Mo95_ppm","Ag107_ppm", "Cd111_ppm", "Sn118_ppm", "Sb121_ppm", "Cs133_ppm", "Ba137_ppm", "W186_ppm", 
  "Tl205_ppm", "Pb207_ppm", "Bi209_ppm", "Th232_ppm", "U238_ppm","Sum_ppm_ICP_MS","TIC", "TN","TC","TS", 
  "TOC", "Microm_001_mean", "Microm_63_mean","Microm_125_mean", "Microm_250_mean", "Microm_500_mean", "Microm_1000_mean")
sample_data <- as.data.frame(unclass(sample_data(ps_ML)))
row.names(sample_data) <- row.names(sample_data(ps_ML))
sample_data2 <- sample_data[, names(sample_data) %in% correlation_cols]
sample_data3 <- sample_data2[, !names(sample_data2) %in% best_correlation_cols]

# adjust naming of input table
subset_list_df <- list(averaged_sample_run10_0_unfiltered.notaxa = sample_data3)

# get response variables
response_variables <- extract_response_variable(
  response_variables = desired_response_vars, phyloseq_object = ps_ML)

# discretize continuous data into two classes
responses <- categorize_response_variable(
  ML_mode = "classification", 
  response_data = response_variables, 
  my_breaks = c(-Inf, 0, Inf),
  class_labels = c("absent", "present"))

# merge the input tables with the response variables
merged_input <- merge_input_response(subset_list_df, responses)

# prepare data for ranger
set.seed(my_seed)
splitted_ranger <- split_data(merged_input, split_ratios)
augmented_ranger <- augment(splitted_ranger, copies, noise_factor)

# prepare data for keras
keras_dummy <- dummify_input_tables(merged_input)
set.seed(my_seed)
splitted_keras <- split_data(keras_dummy, split_ratios)
augmented_keras <- augment(splitted_keras, copies, noise_factor)
scaled_keras <- scaling(augmented_keras)
ready_keras <- inputtables_to_keras(scaled_keras)
# new keras/tensorflow version throws error at first call --> repeat
if(!exists("ready_keras")) {
  ready_keras <- inputtables_to_keras(scaled_keras)
}

# save objects as input for machine learning
saveRDS(augmented_ranger, file = file.path(savePath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
saveRDS(ready_keras, file = file.path(savePath, paste0("input_list_keras_", analysis_ID, ".RDS")))
