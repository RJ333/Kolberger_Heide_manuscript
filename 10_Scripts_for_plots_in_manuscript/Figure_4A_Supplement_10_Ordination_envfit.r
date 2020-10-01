#!/usr/bin/env Rscript

# PCA ordination of unsupervised RF using Top25 genera community composition
# and correlating environmental factors to it
library(ggplot2)
library(dplyr)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# read in supervised ordination data
superv_ordination_table <- read.delim(file.path(resultsPath, 
  "ranger_ross-3_TNT_RF_PCA_genus25_008.tsv"), sep = "\t", header = TRUE)
superv_PCA_variation <- read.delim(file.path(resultsPath, 
  "ranger_ross-3_TNT_RF_PCA_genus25_008_variation.tsv"), sep = "\t", header = TRUE)

# read in unsupervised ordination data
ordination_table <- read.delim(file.path(resultsPath, 
  "ranger_chandler-1_URF_PCA_genus25_008.tsv"), sep = "\t", header = TRUE)
PCA_variation <- read.delim(file.path(resultsPath, 
  "ranger_chandler-1_URF_PCA_genus25_008_variation.tsv"), sep = "\t", header = TRUE)
envfit_results <- read.delim(file.path(resultsPath, 
  "ranger_chandler-1_URF_PCA_genus25_008_envdata.tsv"), sep = "\t", header = TRUE)
envfit_results_sig <- subset(envfit_results, Pvalues < 0.001 & R_squared > 0.3)
envfit_results_sig <- subset(envfit_results_sig, !Variables %in% c("Biological_replicate", 
  "Distance_from_mine", "Weight_loss"))
envfit_results_sig <- envfit_results_sig %>% mutate(Variables = recode(Variables, 
  "Depth_cm" = "Depth", "V51_ppm" = "V", "Co59_ppm" = "Co", "Cu63_ppm" = "Cu", "As75_ppm" = "As", 
  "Mo95_ppm" = "Mo", "Ag107_ppm" = "Ag", "Sb121_ppm" = "Sb", "Cs133_ppm" = "Cs", "Tl205_ppm" = "Tl", 
  "Pb207_ppm" = "Pb", "Bi209_ppm" = "Bi", "U238_ppm" = "U", "TN" = "TN", "TC" = "TC", 
  "TS" = "TS", "TOC" = "TOC", "Microm_001_mean" = ">0.01 µm", "Microm_63_mean" = ">63 µm", 
  "Microm_250_mean" = ">250 µm", "Microm_500_mean" = ">500 µm", "Microm_1000_mean" = ">1000 µm", 
  "Cr52_ppm" = "Cr"))

# plotting settings
extra_scaling_factor <- 2.3
colourblind <- c("#999999", "#CC79A7", "#009E73", "#D55E00", "#F0E442")
point_size <- 2.5
text_size <- 3

# plot Figure 4 A
ordination_envfit <- ggplot(ordination_table, aes(x = X, y = Y)) +
  geom_point(data = subset(ordination_table, Sample_Type == "Surface" & 
      Area %in% c("East_of_restricted", "West_of_restricted")), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(data = subset(ordination_table, Experiment != "Transect" & Area == "Restricted"), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(aes(shape = Target, colour = Area), size = point_size, alpha = 0.7) +
  coord_fixed() +
  geom_segment(data = envfit_results_sig,
               aes(x = 0, xend = Dim1 * extra_scaling_factor * 0.95, y = 0, 
                 yend = Dim2* extra_scaling_factor * 0.95, alpha = R_squared),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_text(data = envfit_results_sig, aes(x = Dim1 * extra_scaling_factor, 
    y = Dim2 * extra_scaling_factor, label = Variables),
            size = text_size, alpha = 0.7, colour = "black") +
  scale_colour_manual(values = colourblind, name = "Areas", 
    breaks = c("Far_west", "West_of_restricted", "Restricted", "Mine_mound", "East_of_restricted"),
      labels = c("Far northwest", "West", "Restricted area", "Mine mound\nin restricted", "East")) +
  scale_alpha_continuous(name = expression(italic(R)^2), breaks = c(0.7, 0.6, 0.5, 0.4)) +
  scale_shape_discrete(name = "TNT", breaks = c("present", "absent"), labels = c("Detected", "Not detected")) +
  guides(color = guide_legend(override.aes = list(size = 2)),
    shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 1),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.75, "lines")) +
  labs(x = paste("PC1 - ", PCA_variation[1,], "%", sep = ""),
    y = paste("PC2 - ", PCA_variation[2,], "%", sep = ""))

# plot Supplement Figure 10, supervised ordination
superv_ordination_envfit <- ggplot(superv_ordination_table, aes(x = X, y = Y)) +
  geom_point(data = subset(superv_ordination_table, Sample_Type == "Surface" & 
      Area %in% c("East_of_restricted", "West_of_restricted")), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(data = subset(superv_ordination_table, Experiment != "Transect" & 
      Area == "Restricted"), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(aes(shape = Target, colour = Area), size = point_size, alpha = 0.7) +
  coord_fixed() +
  scale_colour_manual(values = colourblind, name = "Areas", 
    breaks = c("Far_west", "West_of_restricted", "Restricted", "Mine_mound", "East_of_restricted"),
      labels = c("Far northwest", "West", "Restricted area", "Mine mound\nin restricted", "East")) +
  scale_shape_discrete(name = "TNT", breaks = c("present", "absent"), labels = c("Detected", "Not detected")) +
  guides(color = guide_legend(override.aes = list(size = 2)),
    shape = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 1),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.75, "lines")) +
  labs(x = paste("PC1 - ", superv_PCA_variation[1,], "%", sep = ""),
    y = paste("PC2 - ", superv_PCA_variation[2,], "%", sep = ""))

ggsave(superv_ordination_envfit, file = file.path(plot_path, "supp_Ordination_supervised_genus25.png"),
  device = "png", width = 18, height = 14, dpi = 600, unit = "cm")
