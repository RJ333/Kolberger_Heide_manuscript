#!/usr/bin/env Rscript

# This code displays the robustness of multiple RF and ANN predictions onto
# unsupervised PCA ordination of the Top25 genera community composition. 
library(ggplot2)
library(dplyr)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select the input tables by analysis id and number of trees
analysis_ID <- "ave_run12"
Number_of_trees <- 10000

# get file path for analysis ID and  results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ANN_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ANN", analysis_ID, sep = "_"))
ranger_pred_files <- ranger_files[grepl("stability", ranger_files)]
ranger_pred_files2 <- ranger_pred_files[grepl(Number_of_trees, ranger_pred_files)]
ANN_pred_files <- ANN_files[grepl("stability", ANN_files)]
pred_files <- c(ranger_pred_files2, ANN_pred_files)

# read and combine result frames per ML
predictions <- data.table::rbindlist(lapply(pred_files, read.delim), fill = TRUE)

# read in ordination data
ordination_table <- read.delim(file.path(resultsPath, 
  "ranger_chandler-1_URF_PCA_genus25_008.tsv"), sep = "\t", header = TRUE)
prediction_ordination <- merge(subset(ordination_table, select = -Nodename), 
  predictions, by.x = "Sample", by.y = "Udemm")

# overall prediction accuracy
prediction_ordination$Prediction_rel <- prediction_ordination$Prediction_diff * 100
prediction_ordination$cuts <- cut(prediction_ordination$Prediction_rel, 
  breaks = c(-Inf, 0.5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99.5, Inf), 
  labels = c("<0.5%", "<10%", "<20%", "<30%", "<40%", "<50%", "<60%", "<70%", 
    "<80%", "<90%", "<99.5%", ">99.5%"))
prediction_ordination$stability <- ifelse(prediction_ordination$cuts %in% 
    c("<0.5%", ">99.5%"), "stable", "unstable")
prediction_ordination$wrong <- ifelse(prediction_ordination$cuts == "<0.5%", "correct", "wrong")
prediction_ordination$Y_new <- ifelse(prediction_ordination$ML == "ranger", 
  prediction_ordination$Y, prediction_ordination$Y - 2)

data.table::setDF(predictions)
ranger_predictions <- subset(predictions, ML == "ranger")
ANN_predictions <- subset(prediction_ordination, ML == "ANN", 
  select = c(stability, cuts, Prediction_rel, Udemm, wrong))

### Plot ANN vs RF stability
stability_plot <- ggplot(prediction_ordination, aes(x = X, y = Y_new)) +
 geom_point(data = subset(prediction_ordination, wrong == "correct"), 
   aes(shape = Target, colour = cuts), size = 3, alpha = 0.7) +
 geom_point(data = subset(prediction_ordination, wrong == "wrong"), 
   aes(shape = Target, colour = cuts), size = 3, alpha = 1) +
  scale_colour_manual(drop = FALSE, values = colorRampPalette(
      c("lightblue" ,"green", "yellow", "orange", "red", "brown", "black"))(12),
    name = "Prediction\nwrong") + 
  labs(x = "", y = "") +
  guides(shape = FALSE) +
  annotate(geom = "text", x = - 0.5, y = 1, label = "Random Forest", 
    color = "black", size = 5, fontface = "bold") +
  annotate(geom = "text", x = 2, y = -3, label = "Artificial\nneural networks", 
    color = "black", size = 5, fontface = "bold") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.75, "lines"))

# get "ordination_envfit" from "Figure_4A_Supplement_10..." and combine both plots
ordination_stability <- cowplot::plot_grid(ordination_envfit, stability_plot , 
  labels = c('A', 'B'), label_size = 14, align = "v", axis = "b", ncol = 1, nrow = 2)

# save combined plot
svg(file.path(plot_path, "Ordinations_stability_combined.svg"), 
  width = 18 / 2.54, height = 20 / 2.54)
print(ordination_stability)
dev.off()
ggsave(ordination_stability, file = file.path(plot_path, "Ordinations_stability_combined.png"),
  device = "png", width = 18, height = 20, dpi = 300, unit = "cm")
