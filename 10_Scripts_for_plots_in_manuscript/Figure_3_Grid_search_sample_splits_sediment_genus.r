#!/usr/bin/env Rscript

# Compare Random Forest classifications on validation set for various
# training/test set splits of the input data
library(tidyverse)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select performance metric by column name
performance_metric <- "Balanced_accuracy"
# select the input tables by analysis id
analysis_ID <- "ave_run1_"

# get file path for analysis ID and ML method results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ranger_files <- ranger_files[grepl("processed", ranger_files)]
ranger_files2 <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger_ave_run4"))
ranger_files2 <- ranger_files2[grepl("processed", ranger_files2)]

# read and combine result frames per ML
ranger_sed_results <- data.table::rbindlist(lapply(ranger_files, read.delim), fill = TRUE)
ranger_genus_results <- data.table::rbindlist(lapply(ranger_files2, read.delim), fill = TRUE)
ranger_sed_results_filtered <- ranger_sed_results %>% filter(Number_of_trees == 10000, Mtry_factor == 1)
ranger_genus_results_filtered <- ranger_genus_results %>% filter(Threshold == 0.08, 
  Number_of_trees == 10000, Mtry_factor == 5, Tax_rank == "Genus")

# set taxa factor levels
ranger_sed_results_filtered$Tax_rank <- factor(ranger_sed_results_filtered$Tax_rank, 
  levels = c("notaxa", "Genus"))
ranger_genus_results_filtered$Tax_rank <- factor(ranger_genus_results_filtered$Tax_rank, 
  levels = c("notaxa", "Genus"))

# read and combine result frames per ML
combined_results <- data.table::rbindlist(list(ranger_sed_results_filtered, 
  ranger_genus_results_filtered), fill = TRUE)
combined_results$Nodename <- factor(combined_results$Nodename, 
  levels = c("chandler-1", "joey-2", "ross-3", "monica-4", "rachel-5", "phoebe-6"))
selection <- combined_results %>% 
  filter(Class == "present", step == "training")
selection$analysis_ID2 <- ifelse(selection$analysis_ID == "ave_run1", "Full sediment", 
  "Full community")

# generate statistics on chosen performance metric, requires helper function
ranger_summarized <- selection %>% group_by(ML, Nodename, analysis_ID, mode, step, 
  Number_of_trees, Target, Class, Threshold, Tax_rank, Number_of_samples, Noise_copies, 
  Noise_factor, Mtry_factor, Subset_1, Subset_2, Subset_3) %>% 
  summarize(Mean = mean(.data[[performance_metric]], na.rm = TRUE),
    low05 = confidence_interval(.data[[performance_metric]], 0.05), 
    high95 = confidence_interval(.data[[performance_metric]], 0.95),
    median = median(.data[[performance_metric]], na.rm = TRUE),
    min = min(.data[[performance_metric]], na.rm = TRUE),
    max = max(.data[[performance_metric]],na.rm = TRUE), 
    sd = sd(.data[[performance_metric]], na.rm = TRUE),
    Total_iterations = n_distinct(Cycle))

samples_split <- data.frame(Repetitions = ranger_summarized$Total_iterations, 
  Nodename = ranger_summarized$Nodename, analysis_ID = ranger_summarized$analysis_ID)

# Plot classification per split
splits <- ggplot(selection, aes(x = as.character(as.numeric(Nodename)), y = Balanced_accuracy * 100, colour = analysis_ID)) +
  geom_violin(size = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(data = subset(selection, analysis_ID == "ave_run1"), y = 69, label = "n = 100", 
    size = 3.5, angle = 90, colour = "grey30", position = position_nudge(x = -0.2), hjust = "left") +
  geom_text(data = subset(selection, analysis_ID == "ave_run4"), y = 69, label = "n = 200", 
    size = 3.5, angle = 90, colour = "grey30", position = position_nudge(x = +0.2), hjust = "left") +
  stat_summary(position = position_dodge(width = 0.9), fill = "black", 
    fun = mean, geom = "point", shape = 21, alpha = 1, size = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(68, 85), expand = c(0, 0)) +
    scale_colour_manual(values = c("blue", "red"),
      name = NULL,
      breaks = c("ave_run1", "ave_run4"),
      labels = c("Full sediment", "Full community")) +
  guides(colour = guide_legend(override.aes = list(shape = 19, linetype = 0, size = 2, colour = c("blue", "red")))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.margin=margin(t=0, r=0, b=-0.2, l=0, unit="cm"),
    strip.text = element_text(size = 12, face = "bold")) +
  labs(x = "Data splits", y = "Balanced accuracy [%]")

ggsave(splits, file = file.path(plot_path, "Splits_sediment_community.tiff"),
  device = "tiff", width = 8.5, height = 10, dpi = 300, unit = "cm")
ggsave(splits, file = file.path(plot_path, "Splits_sediment_community.png"),
  device = "png", width = 8.5, height = 10, dpi = 300, unit = "cm")
