#!/usr/bin/env Rscript

# plot false positives classification rate and their metabolite content
library(ggplot2)
library(speedyseq)
library(dplyr)
library(corrplot)

# select project directory, the folder containing the results from classification
# the phyloseq object and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
phyloPath <- file.path(projectPath,"phyloseq_output")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select the input tables by analysis id and number of trees for RF
analysis_ID <- "ave_run12"
Number_of_trees <- 10000
# get file path for analysis ID and ML method results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ANN_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ANN", analysis_ID, sep = "_"))
# select prediction robustness files
ranger_pred_files <- ranger_files[grepl("stability", ranger_files)]
ranger_pred_files2 <- ranger_pred_files[grepl(Number_of_trees, ranger_pred_files)]
ANN_pred_files <- ANN_files[grepl("stability", ANN_files)]
pred_files <- c(ranger_pred_files2, ANN_pred_files)
# read in phyloseq object
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))

# read and combine result frames per ML
predictions <- data.table::rbindlist(lapply(pred_files, read.delim), fill = TRUE)
data.table::setDF(predictions)
sampledata <- as.data.frame(sample_data(ps_ML))
predictions_sample <- merge(predictions, sampledata, by.x = "Udemm", by.y = "row.names", all.x = TRUE)

# correlate concentration with misclassification
misclass_full <- subset(predictions_sample, !is.na(Error_type))
misclass_uxo <- misclass_full[,c("ML", "Udemm", "Error_type", "Prediction_diff",
  "TNT", "ADNT_2", "ADNT_4", "DANT_2.4", "DANT_2.6")]

# select specifically false positives and prepare data for plot
misclass_pos <- subset(misclass_uxo, Error_type == "False_positive")
misclass_pos$Metabolites <- misclass_pos$ADNT_2 + misclass_pos$ADNT_4 + 
  misclass_pos$DANT_2.4 + misclass_pos$DANT_2.6
misclass_pos$Present <- ifelse(misclass_pos$Metabolites > 0, "yes", "no")
misclass_pos$ML_nice <- as.factor(ifelse(misclass_pos$ML == "ranger", "Random Forest", 
  "Artificial Neural Network"))
misclass_pos$ML_nice <- relevel(misclass_pos$ML_nice, "Random Forest")

# plot misclassification of false positives
false_pos_plot <- ggplot(misclass_pos, aes(x = reorder(Udemm, -Prediction_diff), 
  y = Prediction_diff * 100, fill = Present)) +
  geom_point(colour = "white", size = 4, shape = 23, alpha = 1) +
  scale_fill_manual(name = "Metabolites", 
    values = c("black", "red"),
    breaks = c("yes", "no"),
    labels = c("Present", "Absent")) +
  facet_wrap(~ML_nice, scales = "free_x")+
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank()) +
  labs(x = "False positive samples", y = "Prediction error [%]")

ggsave(false_pos_plot, file = file.path(plot_path, "False_positives_error_small.tiff"),
  device = "tiff", width = 18, height = 10, dpi = 300, unit = "cm")
