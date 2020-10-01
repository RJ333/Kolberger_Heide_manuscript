#!/usr/bin/env Rscript

# This script visualizes the ASVs comprising the Top25 genera and compares them
# to those present in control libraries and actual samples
library(speedyseq)
library(ggplot2)
# select project directory,
# the phyloseq objects and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
phyloPath <- file.path(projectPath,"phyloseq_output")
plotPath <- file.path(projectPath, "plots")
# read in the phyloseq objects
ps_full <- readRDS(file.path(phyloPath, "ps_full.RDS"))
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))

# reduce ps_ML to most important 25 genera
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
best_genera <- tax_table(ps_ML)[,"To_Genus"] %in% top25_genera
ps_top25 <- prune_taxa(best_genera, ps_ML)
ps_top25 <- filter_taxa(ps_top25, function (x) {sum(x > 0.08) >= 1}, prune = TRUE)

# subset ps objects to control samples
negative_V4_controls <- subset_samples(ps_full, Library_purpose == "negative_pcr_control" & 
    Primerset == "V4", prune = TRUE)
negative_V4_controls <- filter_taxa(negative_V4_controls, function (x) {sum(x > 0) >= 1}, prune = TRUE)      

positive_V4_controls <- subset_samples(ps_full, Library_purpose == "positive_pcr_control" & 
    Primerset == "V4", prune = TRUE)
positive_V4_controls <- filter_taxa(positive_V4_controls, function (x) {sum(x > 0) >= 1}, prune = TRUE)   
  
blank_V4_controls <- subset_samples(ps_full, Library_purpose == "blank_extraction_control" & 
    Primerset == "V4", prune = TRUE)
blank_V4_controls <- filter_taxa(blank_V4_controls, function (x) {sum(x > 0) >= 1}, prune = TRUE)  

# get vectors of ASVs and compare to ASVs of top 25 genera
in_neg_control <- row.names(tax_table(negative_V4_controls))
in_pos_control <- row.names(tax_table(positive_V4_controls))
in_blank_control <- row.names(tax_table(blank_V4_controls))
in_top25 <- row.names(tax_table(ps_top25))

intersect(in_top25, in_neg_control)
intersect(in_top25, in_pos_control)
intersect(in_top25, in_blank_control)

# plot abundance of control ASVs in the rest of the samples
positive_V4_controls <- subset_samples(ps_full, 
  Library_purpose == "positive_pcr_control" & Primerset == "V4", prune = TRUE)

include <- c("ASV00063", "ASV00074")
intersect(row.names(tax_table(ps_top25)), include)

ps_control_ASV <- prune_taxa(include, ps_full)
ps_control_ASV <- filter_taxa(ps_control_ASV, function (x) {sum(x > 0) >= 1}, prune = TRUE) 
only_controls <- subset_samples(ps_control_ASV, Library_purpose != "sample")
not_controls <- subset_samples(ps_control_ASV, Library_purpose == "sample")
# which genus do they belong to? ASV63 Maribacter (4 reads in positive PCR plate 1) 
# and ASV74 Cobetia (5 reads in negative PCR plate 2), none in blank, up to 3000 reads in samples
tax_table(ps_control_ASV)

# which Samples have the highest abundance of ASVs
plot_bar(not_controls, x = "Library", fill = "OTU") +
  facet_wrap(~ Library_purpose, scales = "free_x")
plot_bar(only_controls, x = "Library", fill = "OTU") +
  facet_wrap(~ Library_purpose, scales = "free_x")
  
ggsave(positive_V4, filename = file.path(plotPath, "positive_V4_cleaned.pdf"), 
  device = "pdf", width = 35, height = 25)
