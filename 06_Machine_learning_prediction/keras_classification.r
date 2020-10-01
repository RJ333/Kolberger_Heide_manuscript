#!/usr/bin/env Rscript
# This script runs keras ANN analyses with the prepared data sets
library(keras)
library(phyloseq2ML)
library(futile.logger)
flog.threshold(INFO)

# set parameters
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

# set settings for ML analysis and environment
ML <- "ANN"
analysis_ID <- "ave_run13"
mode <- "classification"
step <- "training"
nodename <- Sys.info()[["nodename"]]
log_time_mem <- TRUE

# read in input data for machine learning and extract the parameters from the object names
keras_ready <- readRDS(file = file.path(inputPath, paste0("input_list_keras_", analysis_ID, ".RDS")))
parameter_df <- extract_parameters(keras_ready)

# add workflow specific variables
parameter_df$analysis_ID <- analysis_ID
parameter_df$ML <- ML
parameter_df$mode <- mode
parameter_df$step <- step
parameter_df$Nodename <- nodename

# set hyperparameters and combine with input data parameters
hyper_grid <- expand.grid(
  ML_object = names(keras_ready),
  Epochs = 100, 
  Batch_size = 4, 
  Cycle = c(1:100),
  k_fold = 3,                  # set to 1 for prediction
  current_k_fold = 1:3,         # set to 1 for prediction
  Early_callback = "val_loss", #prediction: "accuracy", training: "val_loss"
  Layer1_units = 50,  
  Layer2_units = 20,
  Dropout_layer1 = 0.0,
  Dropout_layer2 = 0.0,
  Dense_activation_function = "relu",
  Output_activation_function = "sigmoid", # sigmoid for binary, softmax for multiclass classification
  Optimizer_function = c("rmsprop", "adam"),
  Loss_function = "binary_crossentropy", # binary_crossentropy for binary, categorical_ for multiclass classification
  Metric = "accuracy",
  Delay = 2)
master_keras <- merge(parameter_df, hyper_grid, by = "ML_object")

# check if settings make sense
if(step == "training") if(unique(hyper_grid$Early_callback) == "val_loss") print("ok") else print("Callback Error!")
if(step == "prediction") if(unique(hyper_grid$Early_callback) == "accuracy") print("ok") else print("Callback Error!")

# order table by current_k_fold 
master_sorted <- master_keras[order(
  master_keras$ML_object,
  master_keras$Cycle, 
  master_keras$current_k_fold
  ), ]
rownames(master_sorted) <- NULL

master_subset <- master_sorted
# use this if the runs should be separated by cycles into subsets due to resource limitations 
#master_subset <- dplyr::filter(master_sorted, Cycle %in% 16:20)  

# run keras analysis using customized "keras_classification_with_log" function
if (log_time_mem) {
  master_subset$results <- purrr::pmap(cbind(master_subset, .row = rownames(master_subset)), 
    keras_classification_with_log, the_list = keras_ready, master_grid = master_subset)
} else {
  master_subset$results <- purrr::pmap(cbind(master_subset, .row = rownames(master_subset)), 
    keras_classification, the_list = keras_ready, master_grid = master_subset)
}
# process results and store
results_df <-  as.data.frame(tidyr::unnest(master_subset, results))
write.table(results_df, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", analysis_ID, "_", part, "_", mode, "_", nodename, "_output_processed_", step, ".tsv")))

part <- "adam_rmsprop" # make sure to update this value if using cycle subsets
