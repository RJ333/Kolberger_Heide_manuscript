#!/usr/bin/env Rscript

# generate pdfs of all important variables and their abundance in supervised (Top25 genera)
# and unsupervised (Top9 sediments) ordination
library(ggplot2)
library(speedyseq)
library(dplyr)

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
ranger_pred_files <- ranger_files[grepl("stability", ranger_files)]
ranger_pred_files2 <- ranger_pred_files[grepl(Number_of_trees, ranger_pred_files)]

# read and combine result frames per ML
predictions <- data.table::rbindlist(lapply(ranger_pred_files2, read.delim), fill = TRUE)
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))
ps_ML2 <- filter_taxa(ps_ML, function (x) {sum(x > 0.08) >= 1}, prune = TRUE) 
ps_ML_genus <- tax_glom(ps_ML2, "To_Genus")

# format labels
sediments <- c(
  "As75_ppm",
  "Microm_63_mean",
  "Microm_250_mean",
  "Microm_500_mean",
  "V51_ppm",
  "Zn66_ppm",
  "TN",
  "Hg_microg_kg",
  "Co59_ppm",
  "Pb207_ppm",
  "Pb206_207_ratio"
  )
sediment_formatted <- c("Arsenic [ppm]", 
  "63 - 125 µm fraction [%]", 
  "250 - 500 µm fraction [%]", 
  "500 - 1000 µm fraction [%]", 
  "Vanadium [ppm]", 
  "Zinc [ppm]", 
  "Total nitrogen [%]", 
  "Mercury [ppb]", 
  "Cobalt [ppm]",
  "Lead [ppm]",
  "Lead isotope ratio [c(206Pb)/c(207Pb)]")
top25_genera <- c(
  "Bacteria_Cyanobacteria_Sericytochromatia_NA_NA_NA",
  "Bacteria_Bacteroidetes_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Maribacter",
  "Bacteria_Proteobacteria_Gammaproteobacteria_Oceanospirillales_Halomonadaceae_Cobetia",  
  "Bacteria_Planctomycetes_Planctomycetacia_Pirellulales_Pirellulaceae_Blastopirellula",  
  "Bacteria_Bacteroidetes_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Maritimimonas",
  "Bacteria_Proteobacteria_Gammaproteobacteria_Oceanospirillales_Nitrincolaceae_Motiliproteus",  
  "Bacteria_Proteobacteria_Gammaproteobacteria_Alteromonadales_Colwelliaceae_NA",
  "Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Thiotrichaceae_Cocleimonas",
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Anaeromicrobium",
  "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfuromonadales_Desulfuromonadaceae_Desulfuromusa",
  "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_MSBL7", 
  "Bacteria_Proteobacteria_Deltaproteobacteria_SAR324_clade(Marine_group_B)_NA_NA",
  "Bacteria_Proteobacteria_Alphaproteobacteria_Caulobacterales_Hyphomonadaceae_NA",  
  "Bacteria_Proteobacteria_Deltaproteobacteria_Bdellovibrionales_Bdellovibrionaceae_OM27_clade",
  "Bacteria_Gemmatimonadetes_Gemmatimonadetes_Gemmatimonadales_Gemmatimonadaceae_NA",
  "Bacteria_Chloroflexi_Anaerolineae_SBR1031_A4b_NA",
  "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_Desulfocapsa",
  "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_Desulfobulbus",
  "Bacteria_Proteobacteria_Alphaproteobacteria_Sphingomonadales_Sphingomonadaceae_Altererythrobacter",
  "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_NA_NA",
  "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodovibrionales_Kiloniellaceae_NA",
  "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfuromonadales_Geobacteraceae_Geopsychrobacter",    
  "Archaea_Euryarchaeota_Methanobacteria_Methanobacteriales_Methanobacteriaceae_Methanobacterium",
  "Bacteria_Planctomycetes_Planctomycetacia_Planctomycetales_Gimesiaceae_NA",
  "Bacteria_TA06_NA_NA_NA_NA")
best_genera <- tax_table(ps_ML_genus)[,"To_Genus"] %in% top25_genera
ps_top25 <- prune_taxa(best_genera, ps_ML_genus)
top25_matrix <- as(otu_table(ps_top25), "matrix")
top25_df <- as.data.frame(top25_matrix)
names(top25_df) <- tax_table(ps_top25)[,"To_Genus"]

# read in ordination data and merge with abundance information of genera and sediment parameters
ordination_table <- read.delim(file.path(resultsPath, # unsupervised ordination
  "ranger_chandler-1_URF_PCA_genus25_008.tsv"), sep = "\t", header = TRUE)
ordination_table_TNT <- read.delim(file.path(resultsPath, # supervised
  "ranger_ross-3_TNT_RF_PCA_genus25_008.tsv"), sep = "\t", header = TRUE)

abu_ordination <- merge(ordination_table_TNT, top25_df, by.x = "Sample", 
  by.y = "row.names")
sedi_ordination <- merge(ordination_table, top25_df, by.x = "Sample", 
  by.y = "row.names")
prediction_abu_ordination <- merge(subset(abu_ordination, select = -Nodename), 
  predictions, by.x = "Sample", by.y = "Udemm")
prediction_sedi_ordination <- merge(subset(sedi_ordination, select = -Nodename), 
  predictions, by.x = "Sample", by.y = "Udemm")

# overall prediction accuracy to percentage
prediction_abu_ordination$Prediction_rel <- prediction_abu_ordination$Prediction_diff * 100
prediction_sedi_ordination$Prediction_rel <- prediction_abu_ordination$Prediction_diff * 100

# plotting concentrations for the most important genera based on
# Top25 genera supervised classification ordination
point_size <- 3
text_size <- 0.7

plotlist <- list()
i <- 0  # plot counter
for (genus in top25_genera) {
  i <- i + 1
  current_letter <- LETTERS[i]
  # automatically decide abundance intervals for coloring
  break_points <- c(0, seq(min(prediction_abu_ordination[[genus]]) - 0.001, 
                      max(prediction_abu_ordination[[genus]]) + 0.001, 
                      length.out = 15))
  nvalues <- sum(table(cut(prediction_abu_ordination[[genus]], break_points, 
    include.lowest = TRUE)) > 0)
  # plot data and add to list
  plotlist[[genus]]  <- ggplot(prediction_abu_ordination, aes(x = X, y = Y)) +
   geom_point(aes(shape = Target, colour = cut(!!sym(genus), !!break_points, 
     include.lowest = TRUE)), size = point_size, alpha = 0.6) +
     coord_fixed() + 
     scale_colour_manual(drop = TRUE, values = colorRampPalette(
        c("grey90", "#0072B2", "#E69F00"))(nvalues),
     name = "Relative\nabundance [%]") + 
    ggtitle(current_letter, subtitle = genus) +
   geom_text(data = subset(prediction_abu_ordination, 
     prediction_abu_ordination[[genus]] > 0), 
      aes(y = Y - 0.02, label = round(!!sym(genus), 2)), colour = "black", 
     size = text_size, alpha = 1) +
   geom_text(data = subset(prediction_abu_ordination, prediction_abu_ordination[[genus]] > 0), 
      aes(y = Y + 0.01, label = round(Prediction_rel, 2)), colour = "black", 
     size = text_size, alpha = 1) + 
    scale_shape_discrete(name = "TNT", breaks = c("present", "absent"), 
      labels = c("Detected", "Not detected")) +
    annotate(geom = "text", x = 2, y = -1.5, label =  "Upper value = Prediction wrong [%]  ", 
      size = text_size + 2)+
    annotate(geom = "text", x = 2, y = -1.7, label = "Lower value = Relative abundance [%]", 
      size = text_size + 2)+
    guides(colour = guide_legend(title = "" ,label = FALSE), 
      shape = guide_legend(override.aes = list(size = 2))) +
     theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle=element_text(size = 8),
      panel.grid = element_blank(),
      legend.title = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.75, "lines")
      ) +
    labs(x =  "", y = "", colour = "")
}
# print list and store in pdf
pdf(file.path(plot_path, "supp_Ordinations_supervised_abundance_genus.pdf"), 
  width = 18 / 2.54, height = 16 / 2.54)
for (plott in plotlist) {
  print(plott)
}
dev.off()

# plotting concentrations for the most important sediment parameters based on
# Top25 genera unsupervised classification ordination
plotlist2 <- list()
j <- 0  # plot counter
for (variable in sediments) {
  j <- j + 1
  current_letter <- LETTERS[j]
  current_sediment <- sediment_formatted[j]
  # automatically decide abundance intervals for coloring
  break_points <- c(0, seq(min(prediction_sedi_ordination[[variable]]) - 0.001, 
                      max(prediction_sedi_ordination[[variable]]) + 0.001, 
                      length.out = 15))
  nvalues <- sum(table(cut(prediction_sedi_ordination[[variable]], break_points, 
    include.lowest = TRUE)) > 0)
  # plot data and add to list
  plotlist2[[variable]]  <- ggplot(prediction_sedi_ordination, aes(x = X, y = Y)) +
   geom_point(aes(shape = Target, colour = cut(!!sym(variable), !!break_points, 
     include.lowest = TRUE)), size = point_size, alpha = 0.7) +
     coord_fixed() + 
     scale_colour_manual(drop = TRUE, values = colorRampPalette(
        c("grey90", "#0072B2", "#E69F00"))(nvalues),
     name = "Concentration") + 
   ggtitle(current_letter, subtitle = current_sediment) +
   geom_text(data = subset(prediction_sedi_ordination, 
     prediction_sedi_ordination[[variable]] > 0), 
      aes(y = Y - 0.02, label = round(!!sym(variable), 2)), colour = "black", 
     size = text_size, alpha = 1) +
   geom_text(data = subset(prediction_sedi_ordination, 
     prediction_sedi_ordination[[variable]] > 0), 
      aes(y = Y + 0.01, label = round(Prediction_rel, 2)), colour = "black", 
     size = text_size, alpha = 1) + 
    scale_shape_discrete(name = "TNT", breaks = c("present", "absent"), 
      labels = c("Detected", "Not detected")) +
    annotate(geom = "text", x = 1.9, y = -1.2, label =  "Upper value = Prediction wrong [%]  ", 
      size = text_size + 2)+
    annotate(geom = "text", x = 1.9, y = -1.4, label = "Lower value = Concentration/Fraction", 
      size = text_size + 2)+
  guides(colour = guide_legend(title = "", label = FALSE), 
    shape = guide_legend(override.aes = list(size = 2))) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle=element_text(size = 8),
      panel.grid = element_blank(),
      legend.title = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.75, "lines")
      ) +
    labs(x =  "", y = "", colour = "")
}
# print list and store in pdf
pdf(file.path(plot_path, "supp_Ordinations_unsupervised_abundance_sediment.pdf"), 
  width = 18 / 2.54, height = 16 / 2.54)
for (plott in plotlist2) {
  print(plott)
}
dev.off()
