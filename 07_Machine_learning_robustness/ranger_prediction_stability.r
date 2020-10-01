#!/usr/bin/env Rscript

# this code allows to record the predictions of several RF models

# set up environment
library(ranger)

# set parameters
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

# set settings for ML analysis and environment
ML <- "ranger"
analysis_ID <- "ave_run12"
step <- "training"
nodename <- Sys.info()[["nodename"]]

# read in input data for machine learning, select genus rank
augmented_input <- readRDS(file = file.path(inputPath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
new_list <- augmented_input[grepl("Genus", names(augmented_input))]

# Random Forest parameters
num.trees <- 1000
repetitions <- c(1:1000)
mtry_factor <- 1

# Run Random Forest
all_predictions <- vector(mode = "list", length = length(new_list))
table_counter <- 0

if (step == "training") {
  for (datatable in new_list) { 
    table_counter <- table_counter + 1
    prediction_results <- vector(mode = "list", length = length(repetitions))
    traintable <- datatable$train_set
    target <- names(traintable)[ncol(traintable)]
    all_vars <- ncol(traintable) - 1
    for_mtry <- ifelse((sqrt(all_vars) * mtry_factor) < all_vars,
      sqrt(all_vars) * mtry_factor, all_vars)
    repetition <- 0
    for (cycle in repetitions) {
      repetition <- repetition + 1
      print(paste("predicting: ", target, "; cycle: ", repetition, "; table: ", table_counter))
      rf_train <- ranger(dependent.variable.name = target, num.trees = num.trees,
        mtry = for_mtry, data = traintable, importance = "none")  
      prediction_results[[paste(target, repetition, table_counter, sep = "_")]] <- rf_train$predictions
    }
    # collect results and calculate average prediction
    prediction_results_df <- as.data.frame(do.call(cbind, prediction_results))
    prediction_results_df[prediction_results_df == 2] <- 0
    prediction_results_df$mean <- rowMeans(prediction_results_df[, c(1:length(repetitions))])
    prediction_results_df <- prediction_results_df[, "mean", drop = FALSE]
    
    all_predictions[[table_counter]] <- cbind(prediction_results_df, traintable[, ncol(traintable), drop = FALSE])
    names(all_predictions)[table_counter] <- paste(names(new_list)[table_counter])
  } 
}

# further process results from RF run
for (i in seq_along(all_predictions)){
  colnames(all_predictions[[i]]) <- c("Mean", "Actual")
}

# calculate prediction accuracies and error types
# 1 is factor level 1 == absent; level 2 (present) was turned to 0
all_predictions2 <- lapply(all_predictions, function(results) {
  results$Actual_num <- ifelse(results$Actual == "absent", 1, 0)
  results$Prediction_diff <- abs(results$Mean - results$Actual_num)
  results$Agree <- ifelse((results$Mean == 0 & results$Actual == "present") | 
    (results$Mean == 1 & results$Actual == "absent"), "Yes", "No")
  results$Error_type <- ifelse(results$Actual == "absent" & results$Mean < 1, "False_positive", "False_negative")
  results$Error_type <- ifelse(results$Agree == "Yes", NA, results$Error_type)
  results
})

# add context data and save file
result_table <- all_predictions2[[1]]
result_table$Udemm <- row.names(result_table)
result_table$Nodename <- nodename
result_table$Number_of_trees <- num.trees
result_table$ML <- ML

write.table(result_table, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", analysis_ID, "_", nodename, "_", step, "_" , num.trees, "_stability.tsv")))
