#!/usr/bin/env Rscript

# this code allows to record the predictions of several ANN models

# set up environment
library(keras)
library(phyloseq2ML)
# set parameters
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "ML/input_tables")
resultPath <- file.path(projectPath, "ML/results")

# set settings for ML analysis and environment
ML <- "ANN"
analysis_ID <- "ave_run12"
step <- "training"
nodename <- Sys.info()[["nodename"]]
my_seed <- readr::parse_number(Sys.info()[["nodename"]])
split_ratios <- 1  # predictions based on training and validation
copies <- 0
noise_factor <- 0.40

# read in input data for machine learning, select genus rank
merged_input <- readRDS(file = file.path(inputPath, paste0("merged_input_list_", analysis_ID, ".RDS")))
new_list <- merged_input[grepl("Genus", names(merged_input))]
keras_dummy <- dummify_input_tables(new_list)
set.seed(my_seed)
splitted_keras <- split_data(keras_dummy, split_ratios)
augmented_keras <- augment(splitted_keras, copies, noise_factor)

# manual scaling of data, exclude dummy variables and response factor column from scaling
train_set <- augmented_keras[[1]]$train_set
response_column <- names(train_set)[ncol(train_set)]
dummy_columns <- names(train_set)[vapply(train_set, phyloseq2ML:::is.dummy, logical(1))]
not_to_scale_columns <- c(dummy_columns, response_column)
# calculate mean and standard deviation from training data
train_mean <- apply(
  train_set[, !names(train_set) %in% not_to_scale_columns], 2, mean)
train_sd <- apply(
  train_set[, !names(train_set) %in% not_to_scale_columns], 2, stats::sd)
# apply mean and SD to train and test data
scaled_train_columns <- scale(
  train_set[, !names(train_set) %in% not_to_scale_columns], 
  center = train_mean, scale = train_sd)
# replace unscaled with scales columns by name
train_set[, colnames(scaled_train_columns)] <- scaled_train_columns

# keras formatting
# separate input data "X" from response variable "y" for keras
X_train <- as.matrix(train_set[, !names(train_set) %in% response_column])
X_train[is.nan(X_train)] <- 0
if (is.factor(train_set[[response_column]])) {
  # reduce level values by 1 to fit python 0-based indexing
  y_train <- keras::to_categorical(as.numeric(train_set[[response_column]]) - 1)
  # add class levels as name attribute to matrix
  attr(y_train, "dimnames")[[2]] <- levels(train_set[[response_column]])
} else {
  y_train <- train_set[[response_column]]
}
training_data <- X_train
training_labels <- y_train
classes <- ncol(training_labels)

# set hyperparameters, total amount of predictions is repetitions * k_fold
Epochs <- 100 
Batch_size <- 4 
repetitions <- c(1:333)
k_fold <- 3                
current_k_folds <- c(1:3)        
Early_callback <- "val_loss"
Layer1_units <- 50  
Layer2_units <- 40
Dropout_layer1 <- 0.0
Dropout_layer2 <- 0.0
Dense_activation_function <- "relu"
Output_activation_function <- "sigmoid" 
Optimizer_function <- c("adam")
Loss_function <- "binary_crossentropy"
Metric <- "accuracy"
Delay <- 2

# lookup to translate between factor levels and class labels
lookup <- stats::setNames(c(colnames(training_labels)), c(0:(classes - 1)))
prediction_results <- vector(mode = "list", length = length(repetitions * k_fold))

# perform machine learning with keras
for (Cycle in repetitions) {
  for (current_k_fold in current_k_folds) {
    indices <- sample(1:nrow(training_data))
    folds <- cut(1:length(indices), breaks = k_fold, labels = FALSE)
    set.seed(Cycle)
    folds2 <- sample(folds)
  
    # split training data into train and validation, by number of folds
    validation_indices <- which(folds2 == current_k_fold, arr.ind = TRUE)    
    validation_data <- training_data[validation_indices, ]
    validation_targets <- training_labels[validation_indices, ]
    partial_train_data <- training_data[-validation_indices, ]
    partial_train_targets <- training_labels[-validation_indices, ]
    
    # build and compile model
    model <- build_the_model(train_data = training_data, classes = classes, Layer1_units = Layer1_units, 
      Layer2_units = Layer2_units, Dropout_layer1 = Dropout_layer1, Dropout_layer2 = Dropout_layer2, 
      Dense_activation_function = Dense_activation_function, Output_activation_function = Output_activation_function, 
      Optimizer_function = Optimizer_function, Loss_function = Loss_function, Metric = Metric)
    
    print(paste("current repetition:", Cycle, "; current kfold:", current_k_fold))
    
    # train model 
    history <- model %>% keras::fit(
      partial_train_data, 
      partial_train_targets,
      epochs = Epochs, 
      batch_size = Batch_size, 
      callbacks = keras::callback_early_stopping(
        monitor = Early_callback,
        patience = Delay,          
        verbose = 0),
      validation_data = list(validation_data, validation_targets),
      verbose = 0)
  
    # predict classes
    val_predictions <- model %>% keras::predict_classes(validation_data)
    
    
    # prepare results
    factor_targets <- phyloseq2ML:::categoric_to_factor(validation_targets)
    predicted <- data.frame(factor_targets, val_predictions)
    
    predicted$Sample <- row.names(validation_data)
    prediction_results[[paste(Cycle, current_k_fold, sep = "_")]] <- predicted
   }
}

# collect results
pred_list <- data.table::rbindlist(prediction_results)
all_predictions <- aggregate(pred_list[, 1:2], list(pred_list$Sample), mean)
row.names(all_predictions) <- all_predictions$Group.1
all_predictions$Group.1 <- NULL
# further process results from ANN predictions
colnames(all_predictions) <- c("Actual_num", "Mean")
all_predictions2 <- all_predictions
all_predictions2$Actual <- lookup[as.character(all_predictions2$Actual_num)]

# calculate prediction accuracies and error types
# for keras, level 1 == TNT present and level 0 == TNT absent (opposite to RF)
all_predictions2$Prediction_diff <- abs(all_predictions2$Mean - all_predictions2$Actual_num)
all_predictions2$Agree <- ifelse((all_predictions2$Mean == 0 & all_predictions2$Actual_num == 0) | 
  (all_predictions2$Mean == 1 & all_predictions2$Actual_num == 1), "Yes", "No")
all_predictions2$Error_type <- ifelse(all_predictions2$Actual_num == 0 & 
  all_predictions2$Mean > 0, "False_positive", "False_negative")
all_predictions2$Error_type <- ifelse(all_predictions2$Agree == "Yes", NA, all_predictions2$Error_type)

# add context data and save file
all_predictions2$Udemm <- row.names(all_predictions2)
all_predictions2$Nodename <- nodename
all_predictions2$ML <- ML

write.table(all_predictions2, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", analysis_ID, "_", nodename, "_", step, "_stability.tsv")))
