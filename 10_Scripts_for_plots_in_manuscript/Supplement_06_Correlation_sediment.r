#!/usr/bin/env Rscript
# This script checks for autocorrelations within the variables of the UDEMM 
# sediment analyses and generates plots from it
library(ggcorrplot)
library(dplyr)
library(phyloseq)

# select project directory, the folder containing the results from classification
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
plot_path <- file.path(projectPath, "paper/figures/raw")
phyloPath <- file.path(projectPath,"phyloseq_output")

# read in phyloseq object and get sample data table
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))
sample_data <- as.data.frame(unclass(sample_data(ps_ML)))
row.names(sample_data) <- row.names(sample_data(ps_ML))
# select variables involved in correlation
# exclude Calcium and Sn due to failed measurement, 
# and Zr because it is too stable for HCl extraction
correlation_cols <- c("DANT_2.6","DANT_2.4","HMX", "RDX", "TNB", "DNB", "TNT", 
  "DNT_2.6", "Tetryl", "DNT_2.4", "ADNT_4", "ADNT_2", "Hg_microg_kg", "Pb206_207_ratio", 
  "P_percent", "Depth_cm", "Fe_percent", "Sc45_ppm","V51_ppm", "Cr52_ppm","Mn55_ppm",
  "Co59_ppm","Ni60_ppm","Cu63_ppm","Zn66_ppm","As75_ppm", "Sr86_ppm","Mo95_ppm",
  "Ag107_ppm", "Cd111_ppm",  "Sb121_ppm", "Cs133_ppm", "Ba137_ppm", "W186_ppm", 
  "Tl205_ppm", "Pb207_ppm", "Bi209_ppm", "Th232_ppm", "U238_ppm","TIC", "TN", "TC", "TS", 
  "TOC", "Microm_001_mean", "Microm_63_mean","Microm_125_mean", "Microm_250_mean", 
  "Microm_500_mean", "Microm_1000_mean")
correlation_data <- sample_data[, names(sample_data) %in% correlation_cols]
correlation_data2 <- correlation_data %>% rename(
  "Core depth" = "Depth_cm",
  "V" = "V51_ppm",
  "Co" = "Co59_ppm",
  "Cu" = "Cu63_ppm",
  "As" = "As75_ppm",
  "Mo" = "Mo95_ppm",
  "Ag" = "Ag107_ppm",
  "Sb" = "Sb121_ppm",
  "Cs" = "Cs133_ppm",
  "Tl" = "Tl205_ppm",
  "Pb" = "Pb207_ppm",
  "Bi" = "Bi209_ppm",
  "U" = "U238_ppm",
  "TN" = "TN",
  "TC" = "TC",
  "TS" = "TS",
  "TOC" = "TOC",
  "TIC" = "TIC",
  "0.01 - 63 µm" = "Microm_001_mean",
  "63 - 125 µm" = "Microm_63_mean",
  "125 - 250 µm" = "Microm_125_mean",
  "250 - 500 µm" = "Microm_250_mean",
  "500 - 1000 µm" = "Microm_500_mean",
  "> 1000 µm" = "Microm_1000_mean",
  "Cr" = "Cr52_ppm" ,
  "2,6-DANT" = "DANT_2.6",
  "2,4-DANT" = "DANT_2.4",
  "HMX" = "HMX",
  "RDX" = "RDX",
  "TNB" = "TNB",
  "DNB" = "DNB",
  "TNT" = "TNT",
  "2,6-DNT" = "DNT_2.6",
  "Tetryl" = "Tetryl",
  "2,4-DNT" = "DNT_2.4",
  "4-ADNT" = "ADNT_4",
  "2-ADNT" = "ADNT_2",
  "Hg" = "Hg_microg_kg",
  "Pb 206/207 ratio" = "Pb206_207_ratio",
  "P" = "P_percent",
  #"Ca" = "Ca_percent",
  "Fe" = "Fe_percent",
  "Sc" = "Sc45_ppm",
  "Mn" = "Mn55_ppm",
  "Ni" = "Ni60_ppm",
  "Zn" = "Zn66_ppm",
  "As" = "As75_ppm",
  "Mo" = "Mo95_ppm",
  "W" = "W186_ppm",
  "Th" = "Th232_ppm",
  #"Zr" = "Zr90_ppm",
  "Cd" = "Cd111_ppm", 
  #"Sn" =  "Sn118_ppm",
  "Sr" = "Sr86_ppm",
  "Ba" = "Ba137_ppm")

# generate correlation information for plots with and without significance
correlation_matrix <- cor(correlation_data2, method = "spearman")
p_values <- cor_pmat(correlation_data2, method = "spearman", exact = FALSE)

# plotting
correlation_plot <- ggcorrplot(correlation_matrix, hc.order = TRUE, 
  outline.col = "white",  type = "full", colors = c("#56B4E9", "white", "#E69F00"), 
  legend.title = "Correlation",p.mat = p_values, sig.level = 0.01) +
  annotate(geom = "text", x = 44, y = 49, size = 5, 
    label = expression(bolditalic("p")*bold(" < 0.01"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)) 

ggsave(correlation_plot, file = file.path(plot_path, "supp_Sediment_correlation.png"),
  device = "png", width = 27, height = 26, dpi = 300, unit = "cm")
