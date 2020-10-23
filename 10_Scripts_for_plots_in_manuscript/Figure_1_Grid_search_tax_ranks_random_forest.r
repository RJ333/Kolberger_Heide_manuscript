#!/usr/bin/env Rscript

# Compare Random Forest classifications on validation set for various
# taxonomic ranks of the input data
library(tidyverse)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select performance metric by column name
performance_metric <- "Balanced_accuracy"
# select the input tables by analysis id
analysis_ID <- "ave_run4"

# get file path for analysis ID and ML method results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ranger_files <- ranger_files[grepl("processed", ranger_files)]
ranger_files2 <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger_ave_run2"))

# read and combine result frames per ML
ranger_results <- data.table::rbindlist(lapply(ranger_files, read.delim), fill = TRUE)
ranger_ASV_results <- data.table::rbindlist(lapply(ranger_files2, read.delim), fill = TRUE)
ranger_ASV_results_filtered <- ranger_ASV_results %>% filter(Threshold == 0.08, 
  Number_of_trees == 10000, Mtry_factor == 5)

# combine ML result frames and rearrange factor level sorting
combined_results <- data.table::rbindlist(list(ranger_results, ranger_ASV_results_filtered), fill = TRUE)
combined_results$Tax_rank <- factor(combined_results$Tax_rank, 
  levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
combined_results$Nodename <- factor(combined_results$Nodename, 
  levels = c("chandler-1", "joey-2", "ross-3", "monica-4", "rachel-5", "phoebe-6"))
selection <- combined_results %>% 
  filter(Class == "present", step == "training")

# generate statistics on chosen performance metric, requires helper function
ranger_summarized <- selection %>% group_by(ML, analysis_ID, mode, step, 
  Number_of_trees, Target, Class, Threshold, Tax_rank, Number_of_samples, 
  Noise_copies, Noise_factor, Mtry_factor, Subset_1, Subset_2, Subset_3) %>% 
  summarize(Mean = mean(.data[[performance_metric]], na.rm = TRUE),
    low05 = confidence_interval(.data[[performance_metric]], 0.05), 
    high95 = confidence_interval(.data[[performance_metric]], 0.95),
    median = median(.data[[performance_metric]], na.rm = TRUE),
    min = min(.data[[performance_metric]], na.rm = TRUE),
    max = max(.data[[performance_metric]],na.rm = TRUE), 
    sd = sd(.data[[performance_metric]], na.rm = TRUE),
    Total_iterations = n_distinct(Cycle))

# get number of samples and variables per taxonomic rank
samples_tax <- data.frame(Repetitions = ranger_summarized$Total_iterations, 
  Tax_rank = ranger_summarized$Tax_rank)
n_variables <- data.frame(Variables = unique(selection$Number_independent_vars), 
  Tax_rank = unique(selection$Tax_rank))

# Plot classification per taxonomic rank
tax_ranks <- ggplot(selection, aes(x = as.factor(Tax_rank), y = Balanced_accuracy * 100)) +
  geom_violin(size = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(data = samples_tax, aes(y = 68.5, label = paste("n =",  Repetitions)), 
    size = 3.5, angle = 90, colour = "grey30", position = position_nudge(x = +0.2), hjust = "left") +
  geom_text(data = n_variables, aes(y = 68.5, label = paste("Taxa =",  Variables)),
    size = 3.5, angle = 90, colour = "grey30", position = position_nudge(x = -0.2), hjust = "left") +
  stat_summary(position = position_dodge(width = 0.9), fill = "black", 
    fun = mean, geom = "point", shape = 21, alpha = 15, size = 2.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(68, 88), expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "none") +
  labs(x = "Taxonomic rank", y = "Balanced accuracy [%]")

ggsave(tax_ranks, file = file.path(plot_path, "Tax_rank_comparison.png"),
  device = "png", width = 8.5, height = 10, dpi = 300, unit = "cm")
ggsave(tax_ranks, file = file.path(plot_path, "Tax_rank_comparison.tiff"),
  device = "tiff", width = 8.5, height = 10, dpi = 300, unit = "cm")
