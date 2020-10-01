#!/usr/bin/env Rscript
# This script runs random forest analyses with the prepared data sets
library(ranger)
library(phyloseq2ML)
library(futile.logger)
flog.threshold(INFO)

# set parameters
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

# set settings for ML analysis and environment
ML <- "ranger"
analysis_ID <- "ave_run11"
mode <- "classification"
step <- "prediction"
nodename <- Sys.info()[["nodename"]]
log_time_mem <- TRUE

# read in input data for machine learning and extract the parameters from the object names
augmented_input <- readRDS(file = file.path(inputPath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
parameter_df <- extract_parameters(augmented_input)

# add workflow specific variables
parameter_df$analysis_ID <- analysis_ID
parameter_df$ML <- ML
parameter_df$mode <- mode
parameter_df$step <- step
parameter_df$Nodename <- nodename

# set hyperparameters and combine with input data parameters
hyper_grid <- expand.grid(
  ML_object = names(augmented_input),
  Number_of_trees = 10000,
  Mtry_factor = 2,
  Importance_mode = c("none"),
  Cycle = 1:500)
master_grid <- merge(parameter_df, hyper_grid, by = "ML_object")

# string arguments needs to be passed as character, not factor level 
master_grid$Target <- as.character(master_grid$Target)

# run ranger analysis using customized "ranger_classification_with_log" function
if (log_time_mem) {
  master_grid$results <- purrr::pmap(cbind(master_grid, .row = rownames(master_grid)), 
      ranger_classification_with_log, the_list = augmented_input, master_grid = master_grid)
} else {
  master_grid$results <- purrr::pmap(cbind(master_grid, .row = rownames(master_grid)), 
      ranger_classification, the_list = augmented_input, master_grid = master_grid)
}
# process results and store
results_df <-  as.data.frame(tidyr::unnest(master_grid, results))
write.table(results_df, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", analysis_ID, "_", mode, "_", nodename, "_output_processed_", step, ".tsv")))
