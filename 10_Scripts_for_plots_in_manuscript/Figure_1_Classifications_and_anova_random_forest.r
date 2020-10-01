#!/usr/bin/env Rscript

# Compare the Random Forest classifications on validation and test set,
# plot the results and calculate p values with ANOVA
library(tidyverse)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# select performance metric by column name
performance_metric <- "Balanced_accuracy"
# select the input tables by analysis id
analysis_ID <- "ave_run"

# get file path for analysis ID and ML method results data frame
ranger_files <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
ranger_files2 <- ranger_files[grepl("processed", ranger_files)]

# read and combine result frames per ML
ranger_results <- data.table::rbindlist(lapply(ranger_files2, read.delim), fill = TRUE)

ranger_run01 <- ranger_results %>% filter(analysis_ID == "ave_run1", Threshold == 0, 
  Mtry_factor == 1, Number_of_trees == 10000)
ranger_run04 <- ranger_results %>% filter(analysis_ID == "ave_run4", Threshold == 0.08, 
  Mtry_factor == 5, Number_of_trees == 10000, Tax_rank %in% c("Genus"))
ranger_run05 <- ranger_results %>% filter(analysis_ID == "ave_run5", Threshold == 0.08, 
  Mtry_factor == 5, Number_of_trees == 10000, Tax_rank %in% c("Genus"))
ranger_run07 <- ranger_results %>% filter(analysis_ID == "ave_run7", Threshold == 0.08, 
  Mtry_factor == 1, Number_of_trees == 5000, Tax_rank %in% c("Genus"))
ranger_run08 <- ranger_results %>% filter(analysis_ID == "ave_run8", Threshold == 0.08, 
  Mtry_factor == 5, Number_of_trees == 10000, Tax_rank %in% c("Genus"))
ranger_run09 <- ranger_results %>% filter(analysis_ID == "ave_run9", Threshold == 0, 
  Mtry_factor == 1, Number_of_trees == 10000)
ranger_run10 <- ranger_results %>% filter(analysis_ID == "ave_run10", Threshold == 0, 
  Mtry_factor == 1, Number_of_trees == 10000)
ranger_run11 <- ranger_results %>% filter(analysis_ID == "ave_run11", Threshold == 0.08, 
  Mtry_factor == 2, Number_of_trees == 10000, Tax_rank %in% c("Genus"))

# combine ML result frames and rearrange factor level sorting
combined_results <- data.table::rbindlist(list(
    ranger_run01, ranger_run04, ranger_run05, ranger_run07, ranger_run08, 
  ranger_run09, ranger_run10, ranger_run11), fill = TRUE)
combined_results$Tax_rank <- factor(combined_results$Tax_rank, 
  levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
combined_results$step <- factor(combined_results$step, 
  levels = c("training", "prediction"))
combined_results$Tax_rank <- fct_explicit_na(combined_results$Tax_rank, "None")
combined_results$Nodename <- factor(combined_results$Nodename, 
  levels = c("chandler-1", "joey-2", "ross-3", "monica-4", "rachel-5", "phoebe-6"))
combined_results$analysis_ID <- factor(combined_results$analysis_ID,
  levels = c("ave_run1", "ave_run9", "ave_run10", "ave_run5", "ave_run11", 
    "ave_run2", "ave_run3", "ave_run4", "ave_run7", "ave_run8"))
# Classes "present" and "absent" are complementary, only choose one
selection <- combined_results %>% filter(Class == "present")

# generate statistics on chosen performance metric, requires helper function "confidence_interval"
all_summarized <- selection %>% 
  group_by(ML, analysis_ID, step, Tax_rank) %>%
  summarize(Mean = mean(.data[[performance_metric]], na.rm = TRUE),
    low05 = confidence_interval(.data[[performance_metric]], 0.05), 
    high95 = confidence_interval(.data[[performance_metric]], 0.95),
    median = median(.data[[performance_metric]], na.rm = TRUE),
    min = min(.data[[performance_metric]], na.rm = TRUE),
    max = max(.data[[performance_metric]],na.rm = TRUE), 
    sd = sd(.data[[performance_metric]], na.rm = TRUE),
    Repetitions = n_distinct(Cycle))

# calculate diff between training and prediction balanced accuracy mean
all_summarized_train <- all_summarized %>% filter(step == "training")
all_summarized_train2 <- all_summarized_train %>% filter(!analysis_ID %in% c("ave_run2", "ave_run3"))
all_summarized_predict <- all_summarized %>% filter(step == "prediction")
all_summarized_predict$Difference <- (all_summarized_train2$Mean - all_summarized_predict$Mean)*100

# plot diff training vs prediction:  positive values: prediction was worse than training
ggplot(all_summarized_predict, aes(x = analysis_ID, y = Difference, colour = ML)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(breaks = c("ave_run1","ave_run2", "ave_run3", "ave_run4", "ave_run5", 
    "ave_run7", "ave_run8", "ave_run9", "ave_run10", "ave_run11"), 
    labels = c("Full sediment" ,"Full community ASV_0.08%", "Full community Genus_0.02%", 
      "Full community Genus_0.08%", "Full combined", "Top25 Genus 0.08%", "Rest Genus 0.08%", 
      "Top9 sediment", "Rest sediment", "Top combined")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1))

# get TNT class ratio averaged per VM for training and prediction
selection$Majority_fraction_class <- ifelse(selection$Positive > selection$Negative, 
  "Positive", "Negative")
selection$Majority_fraction_class_value <- ifelse(selection$Majority_fraction_class == "Positive", 
  selection$Majority_fraction, 100 - selection$Majority_fraction)

selection_train <- subset(selection, step == "training")
table(selection_train$Majority_fraction_class_value)
ratios_train <- c(42.9, 43.8, 44.6, 44.6, 44.6, 48.2)
mean(ratios_train) # 44.8

selection_test <- subset(selection, step == "prediction")
table(selection_test$Majority_fraction_class_value)
ratios_test <- c(36.8, 47.4, 50, 50, 50, 52.6)
mean(ratios_test) # 47.8

#new facet label names
new_steps <- c("Validation set (112 samples)", "Holdout set (38 samples)")
names(new_steps) <- c("training", "prediction")

groups_train <- c("A", "B", "C", "D", "E", "D", "E", "F")
groups_test <- c("G", "H", "I", "J", "J", "K", "L", "M")
 
# plot means of training and test setups
classification_means <- ggplot(all_summarized, aes(x = analysis_ID, y = Mean * 100, colour = analysis_ID)) +
  geom_rect(data = subset(all_summarized, step == "training"), aes(xmin = 0, xmax = 9, 
    ymin = 42.9, ymax = 48.2), alpha = .2, colour = "grey70") +
  geom_rect(data = subset(all_summarized, step == "prediction"), aes(xmin = 0, xmax = 9,
    ymin = 36.8, ymax = 52.6), alpha = .2, colour = "grey70") +
  geom_errorbar(aes(ymin = (Mean - sd) * 100, ymax = (Mean + sd) * 100), width = 0.4, 
    size = 0.8, alpha = 0.4, colour = "grey30") +
  geom_point(alpha = 1, size = 3.5) +
  scale_colour_manual(values = c("blue", "blue", "blue", "blue", "blue", "red", "red", "red")) +
  geom_vline(xintercept = c(3.5, 5.5), colour = "black") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(30, 100), expand = c(0, 0)) +
  geom_text(aes(y = 60, label = paste("n =", Repetitions)), size = 3, angle = 90, colour = "grey30") +
  geom_text(data = subset(all_summarized, step == "training"), aes(y = 92.5), label = groups_train, 
    size = 4, colour = "grey30") +
  geom_text(data = subset(all_summarized, step == "prediction"), aes(y = 92.5), label = groups_test, size = 4, 
    colour = "grey30") +
  geom_label(aes(x = 1.75, y = 45.6), label = "Class distribution", label.size = 0, size = 3, 
    fill = "white", colour = "black") +
  scale_x_discrete(breaks = c("ave_run1","ave_run2", "ave_run3", "ave_run4", "ave_run5", 
    "ave_run7", "ave_run8", "ave_run9", "ave_run10", "ave_run11"), 
    labels = c("Full sediment" ,"Full community ASV_0.08%", "Full community Genus_0.02%", "Full community", 
      "Full combined", "Top25 community", "Non-Top25 community", "Top9 sediment", "Non-Top9 sediment", 
      "Top combined"))+
  labs(x = "Input data",
    y = "Balanced accuracy [%]") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none") +
  facet_wrap(~step, nrow = 1, labeller = labeller(step = new_steps))

svg(file.path(plot_path, "Classification_comparison.svg"), width = 18 / 2.54, height = 14 / 2.54)
print(classification_means)
dev.off()

ggsave(classification_means, file = file.path(plot_path, "Classification_comparison_small.png"),
  device = "png", width = 18, height = 14, dpi = 300, unit = "cm")

# pvalues test between two groups
set1 <- ranger_run07 %>% filter(Class == "present", step == "training")
set2 <- ranger_run11 %>% filter(Class == "present", step == "training")
t.test(set1$Balanced_accuracy, set2$Balanced_accuracy, alternative = "two.sided", var.equal = FALSE)

# anova two way
fulldata <- selection %>% group_by(analysis_ID, step) %>% select(Balanced_accuracy)
# anova one way and one way analysis of means (not assuming equal vairances)
traindata <- fulldata %>% filter(step == "training")
testdata <- fulldata %>% filter(step == "prediction")

res_step.aov <- aov(Balanced_accuracy ~ step, data = fulldata)
res_train.aov <- aov(Balanced_accuracy ~ analysis_ID, data = traindata)
res_test.aov <- aov(Balanced_accuracy ~ analysis_ID, data = testdata)
TukeyHSD(res_step.aov)
TukeyHSD(res_train.aov)  # 4 gegen 5  und 7 gegen 11 nicht sig
TukeyHSD(res_test.aov)   # 5 gegen 11 nicht sig

res_analysis_step.owt <- oneway.test(Balanced_accuracy ~ step, data = fulldata)
res_analysis_train.owt <- oneway.test(Balanced_accuracy ~ analysis_ID, data = traindata)
res_analysis_predict.owt <- oneway.test(Balanced_accuracy ~ analysis_ID, data = testdata)
