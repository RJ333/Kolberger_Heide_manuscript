#!/usr/bin/env Rscript

# analyze variable importance after Altmann et al
library(phyloseq2ML)
library(ranger)

# set static variables
projectPath <- "/data/projects/2019/tnt"
savePath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

ML <- "ranger"
analysis_ID <- "ave_run1"
nodename <- Sys.info()[["nodename"]]

num.trees <- 10000
mtry_factor <- 1
permutations <- 1000

# read in data
augmented_input <- readRDS(file = file.path(savePath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
input_selection <- augmented_input[1]
parameter_df <- extract_parameters(input_selection)
parameter_df$Nodename <- nodename

# ranger parameters
train_data <- input_selection[[1]][["train_set"]]
all_vars <- ncol(train_data) - 1
Target <- names(train_data)[ncol(train_data)]
for_mtry <- ifelse((sqrt(all_vars) * mtry_factor) < all_vars,
                     sqrt(all_vars) * mtry_factor, all_vars)

# run ranger
result_list <- list()
for (iii in 1:50) {
  rf_object <- ranger(dependent.variable.name = Target, importance = "permutation",
    data = train_data, num.trees = num.trees, mtry = for_mtry)
  # extract and select pvalues for variable importance
  altmann_pvalues <- as.data.frame(importance_pvalues(rf_object, method = "altmann", 
    data = train_data, formula = as.formula(paste(Target, " ~ .")), 
    num.permutations = permutations))
  altmann_pvalues$Cycle <- iii
  altmann_pvalues$Variable <- row.names(altmann_pvalues)
  altmann_pvalues$Prediction_error <- rf_object$prediction.error * 100
  print(iii)
  result_list[[iii]] <- altmann_pvalues
}

# rbind processed list items and merge with classification results and taxonomy
altmann_long <- as.data.frame(data.table::rbindlist(result_list))
altmann_long$RF_object <- names(input_selection)
altmann_long_merged <- merge(altmann_long, parameter_df, by.x = "RF_object", by.y = "ML_object")

write.table(altmann_long_merged, sep = "\t", file = file.path(resultPath, 
  paste0("ranger_", analysis_ID, "_", nodename, "_varimpo.tsv")))
