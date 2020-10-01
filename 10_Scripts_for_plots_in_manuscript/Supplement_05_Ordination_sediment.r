#!/usr/bin/env Rscript

# plotting PCA ordination of unsupervised Random Forest using sediment parameters
library(ggplot2)
library(dplyr)
# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
plot_path <- file.path(projectPath, "paper/figures/raw")

# read in ordination data
ordination_table <- read.delim(file.path(resultsPath, 
  "ranger_ross-3_URF_PCA_sediment_full.tsv"), sep = "\t", header = TRUE)
PCA_variation <- read.delim(file.path(resultsPath, 
  "ranger_ross-3_URF_PCA_sediment_full_variation.tsv"), sep = "\t", header = TRUE)

# plotting
colourblind <- c("#999999", "#CC79A7", "#009E73", "#D55E00", "#F0E442")
point_size <- 2.5

ordination_sediment <- ggplot(ordination_table, aes(x = X, y = Y)) +
  geom_point(data = subset(ordination_table, Sample_Type == "Surface" & 
      Area %in% c("East_of_restricted", "West_of_restricted")), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(data = subset(ordination_table, Experiment != "Transect" & Area == "Restricted"), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(aes(shape = Target, colour = Area), size = point_size, alpha = 0.7) +
  coord_fixed() + ## need aspect ratio of 1!
  scale_colour_manual(values = colourblind, name = "Areas", 
    breaks = c("Far_west", "West_of_restricted", "Restricted", "Mine_mound", "East_of_restricted"),
      labels = c("Far northwest", "West", "Restricted area", "Mine mound\nin restricted", "East")) +
  scale_alpha_continuous(name = expression(italic(R)^2), breaks = c(0.7, 0.6, 0.5, 0.4)) +
  scale_shape_discrete(name = "TNT", breaks = c("present", "absent"), 
    labels = c("Detected", "Not detected")) +
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

ggsave(ordination_sediment, file = file.path(plot_path, "supp_Ordination_sediment_full_unsupervised.png"),
  device = "png", width = 18, height = 14, dpi = 600, unit = "cm")
