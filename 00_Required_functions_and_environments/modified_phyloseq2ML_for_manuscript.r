#' This source file contains the manuscript specific versions of functions 
#' which were originally supplied by package phyloseq2ML

#' Run ranger with parameters of data.frame rows.
#'
#' This functions calls ranger using the parameter values in each row of the 
#' provided master_grid, using the data of the list elements. Please have
#' a look at the [ranger doc](https://cran.r-project.org/web/packages/ranger/ranger.pdf)
#' for explanation on the ranger related variables, the arguments are beginning 
#' with "ranger" in the description. Except for `the list`, `master_grid` and `.row`
#' all arguments need to be column names of `master_grid`
#'
#' @param Target char, the response variable
#' @param ML_object factor or char, the name of the corresponding `the_list` item
#' @param Cycle integer, the current repetition
#' @param Number_of_trees ranger, integer, number of trees per forest
#' @param Mtry_factor ranger, factor to multiply default ranger mtry argument
#' @param .row current row of master_grid
#' @param the_list The input tables list
#' @param master_grid the data frame containing all parameter combinations
#' @param step character declaring `training` or `prediction`
#' @param ... further parameters passed on to subfunctions
#'
#' @return a data frame with results and metrics for each row of the master_grid
#'
#' @export
ranger_classification_with_log <- function(master_grid, Target, ML_object, Cycle, 
  Number_of_trees, Mtry_factor, .row, the_list, step, ...) {
  
  if(!all(c("Target", "ML_object", "Cycle", "Number_of_trees", "Mtry_factor") %in% 
      colnames(master_grid))) {
    stop("Ranger parameters do not match column names in master_grid")
  }
  if(is.null(the_list[[ML_object]])) {
    stop("Names in the_list and master_grid do not match")
  }
  if(!is.character(Target)) {
    stop("ranger requires Target as character to work with purr::pmap()")
  }
  if(!is.factor(the_list[[ML_object]][["train_set"]][[Target]])) {
    stop("Response variable is not a factor")
  }
  stopifnot(step == "training" | step == "prediction")
 # print(.row)
  #print(nrow(master_grid))
  state <- paste("Row", .row, "of", nrow(master_grid))
  futile.logger::flog.info(state)
  all_vars <- ncol(the_list[[ML_object]][["train_set"]]) - 1
  # multiply sqrt of variables with Mtry_factor; if greater than available 
  # number of variables, select all variables
  for_mtry <- ifelse((sqrt(all_vars) * Mtry_factor) < all_vars,
    sqrt(all_vars) * Mtry_factor, all_vars)
  n_classes <- length(levels(as.factor(the_list[[ML_object]][["train_set"]][[Target]])))
  # included timing
  train_time <- system.time({
    RF_train <- ranger::ranger(
      dependent.variable.name = Target,  # needs to character, not factor
      data = the_list[[ML_object]][["train_set"]],  # referring to named list item
      num.trees = Number_of_trees,
      mtry = for_mtry,  
      importance = "none")
  })

  if (step == "prediction") {
    prediction_time <- system.time({
      RF_prediction <- stats::predict(object = RF_train, 
        data = the_list[[ML_object]][["test_set"]])
    })  
    confusion_matrix <- table(true = the_list[[ML_object]][["test_set"]][[Target]], 
      predicted = RF_prediction$predictions)
    store_classification_with_log(trained_rf = RF_train, predicted_rf = RF_prediction, 
      confusion_matrix = confusion_matrix, test_set = the_list[[ML_object]][["test_set"]],
      n_classes = n_classes, step = step, train_timings = train_time, prediction_timings = prediction_time)
  } else {  
    store_classification_with_log(trained_rf = RF_train, confusion_matrix = RF_train$confusion.matrix,
      n_classes = n_classes, step = step, train_timings = train_time)
  }
}

#' Store results from ranger classification training and prediction
#'
#' This function extracts information from the ranger objects generated by
#' training or prediction and stores them in a data.frame. It calls the functions 
#' `classification_metrics` and `prediction_accuracy` to generate metrics on
#' classification performance for each class.
#'
#' @param trained_rf the ranger object generated by training with `ranger()`
#' @param predicted_rf the ranger object generated by prediction with `predict()`,
#'   (default: `NULL`)
#' @param confusion_matrix the confusion matrix obtained from 
#'   `trained_rf$confusion.matrix` or generated for `predicted_rf`
#' @param n_classes the number of classes for classification
#' @param step character declaring whether `training` or `prediction` occurs
#' @param ... parameters passed on to `prediction_accuracy`
#' 
#' @return A data frame with one row per ranger run and class
#'
#' @export
store_classification_with_log <- function(trained_rf, predicted_rf = NULL, 
  confusion_matrix, n_classes, step, train_timings, prediction_timings = c(elapsed = 0), ...) {
  
  stopifnot(step == "training" | step == "prediction")
  if(class(trained_rf) != "ranger") {
    stop("trained_rf is not of class ranger")
  }
  if(!is.numeric(n_classes)) {
    stop("n_classes needs to be numeric")
  }
  
  if(class(confusion_matrix) != "table") {
    stop("confusion_matrix is not a table")
  }
  
  results <- data.frame()
  # extract classifications for each class, every class becomes own row
  for (class in 1:n_classes) {
    results[class, "Class"] <- row.names(confusion_matrix)[class] 
    results[class, "True_positive"] <- confusion_matrix[class, class]
    results[class, "False_positive"] <- sum(confusion_matrix[,class]) - 
      confusion_matrix[class,class]
    results[class, "True_negative"] <- sum(confusion_matrix[-class, -class])
    results[class, "False_negative"] <- sum(confusion_matrix[class,]) - 
      confusion_matrix[class,class]
  }  
  # values differing between training and prediction
  if (step == "training") {
    results$Prediction_error <- trained_rf$prediction.error * 100
    results$Number_of_samples <- as.numeric(trained_rf$num.samples)   
  } else {
    results$Prediction_error <- 100 - prediction_accuracy(predicted_rf, ...)   
    results$Number_of_samples <- as.numeric(predicted_rf$num.samples)
  }
  # these values are identical or similar to extract in training and prediction
  results$Variables_sampled <- as.numeric(trained_rf$mtry)
  results$Number_independent_vars <- as.numeric(trained_rf$num.independent.variables)
  results$Tree_type <- trained_rf$treetype
  results$Vars_percent <- as.numeric(results$Variables_sampled / 
    results$Number_independent_vars) * 100
  results <- classification_metrics(results)
  # added for timings and memory
  results$Training_seconds_elapsed <- as.numeric(train_timings[["elapsed"]])
  results$Prediction_seconds_elapsed <- as.numeric(prediction_timings[["elapsed"]])
  results$System_memory_Mb <- as.numeric(try(system(
    "free --mega  | grep ^Mem | tr -s ' ' | cut -d ' ' -f 3", intern = TRUE)))
  #results$System_memory_Mb <- memory.size()  # windows
  results
}

#' Run keras tensorflow classification.
#' 
#' This functions calls keras tensorflow using the parameter values in each row 
#' of the provided master_grid, using the data of the list elements. Please have
#' a look at the keras [fit doc](https://keras.rstudio.com/reference/fit.html)
#' for explanation on the keras related variables, the arguments are beginning 
#' with "keras" in the description. Except for `the list`, `master_grid` and `.row`
#' all arguments need to be column names of `master_grid`
#' 
#' @param Target factor, the response variable
#' @param ML_object factor or char, the name of the corresponding `the_list` item
#' @param Cycle integer, the current repetition
#' @param Epochs keras, integer, how many times should the whole data set be 
#'   passed through the network?
#' @param Batch_size keras, integer, how many samples before updating the weights?
#' @param k_fold integer, the total number of k_folds for cross validation 
#' @param current_k_fold integer, the current k_fold in range 1 : k_fold 
#' @param Early_callback keras, string, a callback metric
#' @param Delay keras, integer, wait for how many epochs before callback happens?
#' @param step character declaring `training` or `prediction`
#' @param the_list The input tables list
#' @param master_grid the data frame containing all parameter combinations
#' @param .row current row of master_grid
#' @param ... additional features passed by pmap call
#'
#' @return a compiled keras sequential model with two hidden layers
#'
#' @export
keras_classification_with_log <- function(Target, ML_object, Cycle, Epochs, Batch_size, k_fold, 
  current_k_fold, Early_callback, Delay, step, the_list, master_grid, .row, ...) {
  
  if(!all(c("Target", "ML_object", "Cycle", "Epochs", "Batch_size", "k_fold", 
    "current_k_fold", "Early_callback", "Delay", "step") %in% colnames(master_grid))) {
    stop("Keras parameters do not match column names in master_grid")
  }
  if(is.null(the_list[[ML_object]])) {
    stop("Names of items in the_list and ML_object in master_grid do not match")
  }
  if(!exists(c("trainset_labels", "trainset_data", "testset_labels", 
    "testset_data"), where = the_list[[1]])) {
    stop("Item in the_list does not have all required elements:
      trainset_labels, trainset_data, testset_labels, testset_data")
  }
  stopifnot(step == "training" | step == "prediction")

  state <- paste("Row", .row, "of", nrow(master_grid))
  futile.logger::flog.info(state)
  community_table <- the_list[[ML_object]]
  training_data <- community_table[["trainset_data"]]
  training_labels <- community_table[["trainset_labels"]]
  classes <- ncol(training_labels)
  if(classes < 2) {
    stop("Less then 2 classes found, response variable setup seems incorrect")
  }
  # lookup to translate between factor levels and class labels
  lookup <- stats::setNames(c(colnames(training_labels)), c(0:(classes - 1)))
  
  if (step == "prediction" & (k_fold != 1 | current_k_fold != 1)) {
    stop("k_fold and current_k_fold need to be 1 for prediction")
  } else if (step == "training") {
    indices <- sample(1:nrow(training_data))
    folds <- cut(1:length(indices), breaks = k_fold, labels = FALSE)
    set.seed(Cycle)
    folds2 <- sample(folds)
  }
  
  if (step == "training") {

    kfold_msg <- paste("k_fold", current_k_fold, "of", k_fold)
    futile.logger::flog.info(kfold_msg)
    # split training data into train and validation, by number of folds
    validation_indices <- which(folds2 == current_k_fold, arr.ind = TRUE)    
    validation_data <- training_data[validation_indices, ]
    validation_targets <- training_labels[validation_indices, ]
    partial_train_data <- training_data[-validation_indices, ]
    partial_train_targets <- training_labels[-validation_indices, ]
    train_time <- system.time({
      # build and compile model
      model <- build_the_model(train_data = training_data, classes = classes, ...)
      
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
    })
  } else if (step == "prediction") {
    validation_data <- community_table[["testset_data"]]
    validation_targets <- community_table[["testset_labels"]]
    partial_train_data <- training_data
    partial_train_targets <- training_labels
    train_time <- system.time({
      # build and compile model
      model <- build_the_model(train_data = training_data, classes = classes, ...)
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
        test_split = 0.0,
        verbose = 0)
    })
  }
  
  # predict classes
  prediction_time <- system.time({
    val_predictions <- model %>% keras::predict_classes(validation_data)
  })
  
  # prepare results
  factor_targets <- phyloseq2ML:::categoric_to_factor(validation_targets)
  predicted <- data.frame(factor_targets, val_predictions)
  predicted_labels <- data.frame(lapply(predicted, function(i) 
    lookup[as.character(i)]))
  if (nrow(predicted_labels) != nrow(validation_data)) {
    stop("Length of predictions and data to be predicted differs")
  }
  # provide all classes as factor levels, otherwise confusion matrix breaks if
  # a class is not predicted or present at all
  predicted_labels$val_predictions <- factor(predicted_labels$val_predictions, 
    levels = colnames(training_labels))
  predicted_labels$factor_targets <- factor(predicted_labels$factor_targets, 
    levels = colnames(training_labels))
  # calculate confusion matrix
  confusion_matrix <- table(
    true = predicted_labels$factor_targets,
    predicted = predicted_labels$val_predictions)

  # return results data.frame
  store_classification_results_with_log(hist = history,
    prediction_table = predicted_labels, confusion_matrix = confusion_matrix, 
    train_data = partial_train_data, n_classes = classes, 
    train_timings = train_time, prediction_timings = prediction_time)
}

#' Store results from keras tf classification training and prediction
#'
#' This function extracts information from the keras model generated by training
#' or prediction and stores them in a data.frame. By calling `classification_metrics` 
#' various metrics for classification performance are calculated for each class.
#'
#' @param hist the keras history object
#' @param prediction_table the data.frame comparing predictions and true values
#' @param n_classes the number of classes for classification
#' @param confusion_matrix the confusion matrix generated from `prediction_table`
#' @param train_data the training set data.frame
#' 
#' @return A data frame with one row per keras run and class
#'
#' @export
store_classification_results_with_log <- function(hist, prediction_table, n_classes,
  confusion_matrix, train_data, train_timings, prediction_timings, ...) {
  
  if(!is.data.frame(prediction_table)) {
    stop("prediction table is not a data frame")
  } else if(nrow(prediction_table) == 0) {
    stop("prediction table is empty")
  }
  
  results <- data.frame()
  
  # extract classifications for each class, every class becomes own row
  if(nrow(confusion_matrix) == ncol(confusion_matrix)) {
    for (class in 1:n_classes) {
      results[class, "Class"] <- row.names(confusion_matrix)[class] 
      results[class, "True_positive"] <- confusion_matrix[class, class]
      results[class, "False_positive"] <- sum(confusion_matrix[, class]) - 
        confusion_matrix[class, class]
      results[class, "True_negative"] <- sum(confusion_matrix[-class, -class])
      results[class, "False_negative"] <- sum(confusion_matrix[class, ]) - 
        confusion_matrix[class, class]
    }
  } else {
    print("Confusion matrix with odd dimensions")
  }
  
  results$Number_of_samples_train <- nrow(train_data)
  results$Number_of_samples_validate <- nrow(prediction_table)
  results$Number_independent_vars <- ncol(train_data)
  if(nrow(confusion_matrix) == ncol(confusion_matrix)) {
    results <- classification_metrics(results, results$Number_of_samples_validate)
  }
  # added for timings and memory
  results$Training_seconds_elapsed <- as.numeric(train_timings[["elapsed"]])
  results$Prediction_seconds_elapsed <- as.numeric(prediction_timings[["elapsed"]])
  results$System_memory_Mb <- as.numeric(try(system(
    "free --mega  | grep ^Mem | tr -s ' ' | cut -d ' ' -f 3", intern = TRUE)))
  results
}

#' #' Remove ASV/OTU/etc from a phyloseq object based on relative abundance.
#'
#' This version of the filter_taxa function can be used with lists and 
#' lapply or for loops.
#'
#' @param phyloseq_object the object to be filtered
#' @param threshold this relative abundance in percentage have to appear for an ASV, 
#'   than all ASV reads are kept
#' @param num_samples the number of times (samples) the threshold has to be met
#'
#' @return The subsetted phyloseq object
#'
#' @export
filter_subsets_abs <- function(phyloseq_object, threshold, num_samples = 1) {
  if (threshold <= 0) {
    stop("threshold is <= 0 %, filter would have no effect, 
      stopped function")
   # } else if (threshold >= 100) {
    # stop("threshold is >= 100 % relative abundance, stopped function")
  }
  if (num_samples <= 0)
    stop("number of samples is not larger 0, filter would have no effect, 
      stopped function")
  
  phyloseq_subset2 <- phyloseq::filter_taxa(phyloseq_object, function (x) {
    sum(x > threshold) >= num_samples }, prune = TRUE)
  return(phyloseq_subset2)
}

#' Create subsets of phyloseq objects based on ASV/OTU counts and taxa level.
#'
#' This subsetting function allows to create multiple combinations of subsets 
#' from a list of phyloseq objects. The lowest taxonomic level 
#' (usually ASV or OTU) is always included, further taxonomic levels can be 
#' specified as described below. The call to `phyloseq::tax_glom` can be slow, 
#' therefore package `speedyseq` is used, if installed.
#'
#' @param subset_list a list of phyloseq objects
#' @param thresholds an integer vector specifying the input
#'   to `filter_taxa(sum(x > threshold) >= num_samples)`, will be formatted with
#'   with digits after . later string splits
#' @param tax_ranks specifying the tax ranks to agglomerate in the form 
#'   of `setNames(c("To_genus", To_family), c("Genus", "Family"))`. 
#'   Here, "To_genus" is the corresponding taxonomic level in tax_table() and 
#'   "Genus" is appended to the name the agglomerated data.frame in the
#'   results list for later distinction. Check taxa rank using 
#'   `colnames(tax_table(TNT_communities))`
#' @param taxa_prefix The leading name of your taxa, e.g. `ASV` or `OTU`, 
#'   must not contain an underscore or white space
#' @param ... further argument passed on to filter_subsets()
#'
#' @return A list of subsetted community tables for each combination of 
#'   phyloseq_subset, thresholds and tax_ranks (+ ASV/OTU)
#' @export
create_community_table_subsets_abs <- function(subset_list, thresholds, 
  tax_ranks = NULL, taxa_prefix, ...) {
  if (!is.list(subset_list))
    stop("Input needs to be a list")
  if (length(thresholds) < 1)
    stop("No count thresholds provided for subsetting")
  if (!is.character(taxa_prefix) | !length(taxa_prefix) == 1)
    stop("Please provide a single character string as name")
  if(grepl("_", taxa_prefix, fixed = TRUE))
     stop("taxa_prefix needs to be a string without underscores")
  
  subset_list_filtered <- list()
  filter_counter <- 0
  for (phyloseq_subset in subset_list) {
    filter_counter <- filter_counter + 1 
    for (threshold in thresholds) { 
      current_name <- paste(names(subset_list)[filter_counter], 
        threshold, "filtered", sep = "_")
      message <- paste("Generating subset:", current_name)
      futile.logger::flog.info(message)
      subset_list_filtered[[current_name]] <- filter_subsets_abs(phyloseq_subset, threshold, ...)
    }
  }
  if (is.null(tax_ranks)) {
    futile.logger::flog.info("No taxonomic levels for agglomeration specified")
    names(subset_list_filtered) <- paste0(names(subset_list_filtered), ".", taxa_prefix)
    subset_list_comb <- subset_list_filtered
  } else {
    if(!requireNamespace("speedyseq", quietly = TRUE)) {
      futile.logger::flog.info("Applying taxonomic agglomeration, speed could be improved by installing speedyseq")
      subset_list_filtered_tax <- unlist(lapply(subset_list_filtered, function(xx) {
        lapply(tax_ranks, function(yy) {
          phyloseq::tax_glom(xx, taxrank = yy)
        }
        )}
      ))
    } else {
      futile.logger::flog.info("Applying taxonomic agglomeration")
      subset_list_filtered_tax <- unlist(lapply(subset_list_filtered, function(xx) {
        lapply(tax_ranks, function(yy) {
          speedyseq::tax_glom(xx, taxrank = yy)
        }
        )}
      ))
    }  
    names(subset_list_filtered) <- paste0(names(subset_list_filtered),  ".", taxa_prefix)
    subset_list_comb <- c(subset_list_filtered, subset_list_filtered_tax)
  }
  subset_list_comb
}
