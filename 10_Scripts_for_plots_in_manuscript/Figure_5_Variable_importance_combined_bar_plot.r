#!/usr/bin/env Rscript

# plot variable importance of Top25 genus and Top9 sediment
library(speedyseq)
library(tidyverse)
library(ggtext)
library(glue)

# select project directory, the folder containing the results from classification
# the phyloseq object and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
resultsPath <- file.path(projectPath, "ML_results")
phyloPath <- file.path(projectPath,"phyloseq_output")
plot_path <- file.path(projectPath, "paper/figures/raw")
# select the input tables by analysis id
analysis_ID <- "ave_run4"

# read in phyloseq object and retrieve sample data from it
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))
sample_data <- as.data.frame(sample_data(ps_ML))

# get file path for analysis ID and select importance tables
ranger_files_sed <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger","ave_run1", sep = "_"))
ranger_files_genus <- list.files(path = resultsPath, full.names = TRUE, 
  pattern = paste("ranger", analysis_ID, sep = "_"))
impo_files_sed <- ranger_files_sed[grepl("varimpo", ranger_files_sed)]
impo_files_genus <- ranger_files_genus[grepl("impo", ranger_files_genus)]

# sediment variables: read and combine result frames per ML
results_sed <- data.table::rbindlist(lapply(impo_files_sed, read.delim), fill = TRUE)
sediment_per_VM <- results_sed %>% group_by(Nodename, Variable) %>% 
  summarise(importance = mean(importance), pvalue = mean(pvalue))
sediment <- sediment_per_VM %>% group_by(Variable) %>% 
  summarize(importance = mean(importance), pvalue = mean(pvalue)) %>% arrange(-importance)
sediment_best <- subset(sediment, importance > 0.009)
sediment_best$Input <- "Full sediment"
mean_prediction_sed <- results_sed %>% summarise(Prediction_error_mean = mean(Prediction_error))

# genus variables: read and combine result frames per ML
results_genus <- data.table::rbindlist(lapply(impo_files_genus, read.delim), fill = TRUE) %>% 
  filter(Tax_rank %in% c("Genus"))
genus_per_VM <- results_genus %>% group_by(Nodename, Variable) %>% 
  summarise(importance = mean(importance), pvalue = mean(pvalue))
genus <- genus_per_VM %>% group_by(Variable) %>% 
  summarise(importance = mean(importance), pvalue = mean(pvalue))
genus_best <- genus %>% filter(pvalue < 0.05, importance > 0.25)
genus_best$Input <- "Full community"
mean_prediction <- results_genus %>% summarise(Prediction_error_mean = mean(Prediction_error))

# subset ps object based on ASV, includes To_Genus information
ASVs <- unique(genus_best$Variable)
ASV_vector <- row.names(tax_table(ps_ML)) %in% ASVs 
ps_best <- prune_taxa(ASV_vector, ps_ML)  # tax_table(ps_best)[, "To_Genus"]
ASV_taxonomy <- tax_table(ps_best)[,c(6,8,10,11,12)]

# create new column "Variable" which either contains taxonomic or sediment variable name
genus_best_tax <- merge(ASV_taxonomy, genus_best, by.x = "row.names", by.y = "Variable", all.y = TRUE)
genus_best_tax$Variable <- ifelse(is.na(genus_best_tax$To_Genus) & is.na(genus_best_tax$Genus), 
  as.character(genus_best_tax$Row.names), as.character(genus_best_tax$To_Genus))
genus_best_tax2 <- subset(genus_best_tax, pvalue < 0.05 & importance > 0.15)

combined <- data.table::rbindlist(list(genus_best_tax2, sediment_best), fill = TRUE)

# format labels
original_names <- c("Bacteria_Cyanobacteria_Sericytochromatia_NA_NA_NA",                                               
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
  "Bacteria_TA06_NA_NA_NA_NA",
  "As75_ppm",
  "Microm_63_mean",
  "Microm_250_mean",
  "Microm_500_mean",
  "V51_ppm",
  "Zn66_ppm",
  "TN",
  "Hg_microg_kg",
  "Co59_ppm") 

modified <-  c(
  glue("  unclassified *Sericytochromatia* "),                                               
  glue("  *Maribacter* "),               
  glue("  *Cobetia* "),           
  glue("  *Blastopirellula* "),             
  glue("  *Maritimimonas* "),            
  glue("  *Motiliproteus* "),     
  glue("  unclassified *Colwelliaceae* "),                   
  glue("  *Cocleimonas* "),          
  glue("  *Anaeromicrobium* "),                     
  glue("  *Desulfuromusa* "),
  glue("  *Desulfobulbaceae* MSBL7 "),            
  glue("  *Deltaproteobacteria* SAR324/MGB "),                  
  glue("  unclassified *Hyphomonadaceae* "),                  
  glue("  *Bdellovibrionaceae* OM27 clade "),     
  glue("  unclassified *Gemmatimonadacea* "),                
  glue("  *Anaerolineae* SBR1031 A4b "),                                                
  glue("  *Desulfocapsa* "),    
  glue("  *Desulfobulbus* "),   
  glue("  *Altererythrobacter* "),
  glue("  unclassified *Bacteroidales* "),                                         
  glue("  unclassified *Kiloniellaceae* "),                 
  glue("  *Geopsychrobacter* "),
  glue("  *Methanobacterium* "),
  glue("  unclassified *Gimesiaceae* "),                        
  glue("  *Bacteria* TA06 "),
  glue("  Arsenic "),
  glue("  63 - 125 µm fraction "),
  glue("  250 - 500 µm fraction "),
  glue("  500 - 1000 µm fraction "),
  glue("  Vanadium "),
  glue("  Zinc "),
  glue("  Total nitrogen "),
  glue("  Mercury "),
  glue("  Cobalt "))

# plot sorted and averaged importanced per variable
combined_plot <- ggplot(combined, aes(x = reorder(Variable, importance), y = importance)) +
  geom_bar(aes(fill = pvalue), stat = "identity") +
  scale_fill_gradient(name = expression(bolditalic(" p")*bold(" value")), 
    low = "lightblue", high = "black", breaks = c(0.000, 0.005, 0.01)) +

  coord_flip() +
  scale_x_discrete(breaks = original_names, labels = modified) +
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_markdown(),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = c(0.85, 0.80),
    legend.margin=margin(t = 0.2, r = 0.2, b = 0.3, l = 0.2, unit = "cm"),
    legend.background = element_rect(fill = "grey90",
                                  size = 0.5, linetype = "solid", 
                                  colour = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = NULL,  y = "Variable importance") +
  facet_wrap(~Input, scales = "free_x", ncol = 2)

ggsave(combined_plot, file = file.path(plot_path, "Variable_importance_combined.tiff"),
  device = "tiff", width = 18, height = 16, dpi = 300, unit = "cm")
ggsave(combined_plot, file = file.path(plot_path, "Variable_importance_combined.png"),
  device = "png", width = 18, height = 16, dpi = 300, unit = "cm")
