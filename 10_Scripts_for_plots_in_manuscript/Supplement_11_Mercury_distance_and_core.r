#!/usr/bin/env Rscript

# plot mercury concentrations and distribution
library(speedyseq)
library(ggplot2)
library(dplyr)

# select project directory,
# the phyloseq object and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
phyloPath <- file.path(projectPath,"phyloseq_output")
plot_path <- file.path(projectPath, "paper/figures/raw")
# read in phyloseq object
ps_final <- readRDS(file.path(phyloPath, "ps_final.RDS"))
# get sample data
sampledata_full <- sample_data(ps_final)
sampledata167 <- unique(sampledata_full[, -c(11,21:30)])
sampledata2 <- subset(sampledata167, Distance_from_mine < 10)
# adjust MUC data for plotting
cores <- subset(sampledata167, Sample_Type == "Core")
levels(cores$Station_ID) <- c("Profile 1", "Profile 2", "Profile 1", "Profile 2")
levels(cores$Area) <- c("East", "West")
cores$Area <- relevel(cores$Area, "West")

# distribution of Hg and Pb depending on distance from a mine
sampledata2 %>% group_by(Distance_from_mine) %>% add_tally() %>% summarize(
  mean = mean(Hg_microg_kg), median = median(Hg_microg_kg), n = mean(n))
sampledata2 %>% group_by(Distance_from_mine) %>% add_tally() %>% summarize(
  mean = mean(Pb207_ppm), median = median(Pb207_ppm), n = mean(n))

# Hg concentrations all sediments
distance_plot <- ggplot(sampledata2, aes(x = Distance_from_mine, y = Hg_microg_kg)) +
  stat_summary(colour = "black", fun = median, geom = "point", size = 3.5, alpha = 1) +
  geom_jitter(shape = 21, colour = "black", fill = "white", alpha = 0.7, width = 0.1, height = 0.00) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.75, "lines")) +
  labs(x = "Distance from mine [m]",  y = "Mercury [ppb]")

ggsave(distance_plot, file = file.path(plot_path, "supp_mercury_distance.png"),
  device = "png", width = 16, height = 14, dpi = 300, unit = "cm")

# Hg concentrations in MUC samples
mercury_core <- ggplot(cores, aes(x = Hg_microg_kg, y = Depth_cm, colour = Station_ID)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_y_reverse(breaks = scales::pretty_breaks(n = 10)) +
  facet_wrap(~ Area, scales = "free_x", dir = "h") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")) +
  labs(x = "Mercury [ppb]",  y = "Depth [cm]")

ggsave(mercury_core, file = file.path(plot_path, "supp_mercury_core.png"),
  device = "png", width = 16, height = 14, dpi = 300, unit = "cm")
