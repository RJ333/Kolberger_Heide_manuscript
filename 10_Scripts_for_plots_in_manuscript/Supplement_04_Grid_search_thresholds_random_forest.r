#!/usr/bin/env Rscript

# Compare Random Forest classifications on validation set for various
# relative abundance thresholds of the input data

# analyze the results' data frame
library(tidyverse)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select performance metric by column name
performance_metric <- "Balanced_accuracy"
# select the input tables by analysis id
analysis_ID <- "ave_run2"

# get file path for analysis ID and ML method results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ranger_results <- data.table::rbindlist(lapply(ranger_files, read.delim), fill = TRUE)

# combine ML result frames and rearrange factor level sorting
combined_results <- data.table::rbindlist(list(ranger_results), fill = TRUE)
combined_results$Tax_rank <- factor(combined_results$Tax_rank, 
  levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
combined_results$Nodename <- factor(combined_results$Nodename, 
  levels = c("chandler-1", "joey-2", "ross-3", "monica-4", "rachel-5", "phoebe-6"))
selection <- combined_results %>% 
  filter(Class == "present", step == "training", Mtry_factor == "5", Number_of_trees == "10000")

# generate statistics on chosen performance metric, requires helper function
ranger_summarized <- selection %>% group_by(ML, analysis_ID, mode, step, Number_of_trees, 
  Target, Class, Threshold, Tax_rank, Number_of_samples, Noise_copies, 
  Noise_factor, Mtry_factor, Subset_1, Subset_2, Subset_3) %>% 
  summarize(Mean = mean(.data[[performance_metric]], na.rm = TRUE),
    low05 = confidence_interval(.data[[performance_metric]], 0.05), 
    high95 = confidence_interval(.data[[performance_metric]], 0.95),
    median = median(.data[[performance_metric]], na.rm = TRUE),
    min = min(.data[[performance_metric]], na.rm = TRUE),
    max = max(.data[[performance_metric]],na.rm = TRUE), 
    sd = sd(.data[[performance_metric]], na.rm = TRUE),
    Total_iterations = n_distinct(Cycle))

# Plot classification per threshold
thresholds <- ggplot(selection, aes(x = as.factor(Threshold), y = Balanced_accuracy * 100)) +
  geom_violin(size = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(y = 72.5, label = "n = 50", size = 4, angle = 90, colour = "grey30") +
  stat_summary(position = position_dodge(width = 0.9), fill = "black", 
    fun = mean, geom = "point", shape = 21, alpha = 1, size = 3.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(70, 90), expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "none") +
  labs(x = "Relative abundance\nthreshold [%]", y = "Balanced accuracy [%]")

ggsave(thresholds, file = file.path(plot_path, "supp_Abundance_filter_comparison_small.png"),
  device = "png", width = 18, height = 14, dpi = 300, unit = "cm")
