#!/usr/bin/env Rscript

# analyze variable importance after Janitza et al
library(phyloseq2ML)
library(ranger)

# static variables
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

ML <- "ranger"
analysis_ID <- "ave_run5"
mode <- "binary_class"
nodename <- Sys.info()[["nodename"]]
tax_rank <- "Genus"
num.trees <- 10000
mtry_factor <- 5

# read files
augmented_input_full <- readRDS(file = file.path(inputPath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
augmented_input <- augmented_input_full[grepl(tax_rank, names(augmented_input_full))]
parameter_df <- extract_parameters(augmented_input)

# ranger parameters
train_data <- augmented_input[[1]][["train_set"]]
all_vars <- ncol(train_data) - 1
target <- names(train_data)[ncol(train_data)]
for_mtry <- ifelse((sqrt(all_vars) * mtry_factor) < all_vars,
                   sqrt(all_vars) * mtry_factor, all_vars)
# run ranger
result_list <- vector(mode = "list", length = 100)
for (iii in 1:100) {
  rf_object <- ranger(dependent.variable.name = target, importance = "impurity_corrected",
                    data = train_data, num.trees = num.trees, mtry = for_mtry)
  
  # extract and select pvalues for variable importance
  janitza_pvalues <- as.data.frame(importance_pvalues(rf_object, method = "janitza"))
  janitza_pvalues$Cycle <- iii
  janitza_pvalues$Variable <- row.names(janitza_pvalues)
  janitza_pvalues$Prediction_error <- rf_object$prediction.error * 100
  print(iii)
  result_list[[iii]] <- janitza_pvalues
}

janitza_long <- as.data.frame(data.table::rbindlist(result_list))
janitza_long <- subset(janitza_long, pvalue < 0.05, importance > 0.01)
janitza_long$Nodename <- nodename
janitza_long$RF_object <- names(augmented_input)
janitza_long_merged <- merge(janitza_long, parameter_df, by.x = "RF_object", by.y = "ML_object")

write.table(janitza_long_merged, sep = "\t", file = file.path(resultPath, 
  paste0("ranger_", analysis_ID, "_", tax_rank, "_", nodename, "_communityimpo.tsv")))
