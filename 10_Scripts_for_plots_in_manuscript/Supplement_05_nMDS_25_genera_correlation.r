#!/usr/bin/env Rscript

# plot nMDS ordination with correlating environmental factors for 25 genera
library(speedyseq)
library(dplyr)
library(ggplot2)
library(vegan)

# select project directory,
# the phyloseq object and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
phyloPath <- file.path(projectPath,"phyloseq_output")
plot_path <- file.path(projectPath, "paper/figures/raw")

# read in phyloseq object and select 25 most important genera
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))
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
ps_ML2 <- filter_taxa(ps_ML, function (x) {sum(x > 0.08) >= 1}, prune = TRUE)  
best_genera <- tax_table(ps_ML2)[,"To_Genus"] %in% top25_genera
ps_top25 <- prune_taxa(best_genera, ps_ML2)
ps_top25_genus <- tax_glom(ps_top25, "To_Genus")

# get data from phyloseq object
V4_mat <- as(otu_table(ps_top25_genus), "matrix")
V4_otu <- as.data.frame(V4_mat)
V4_sample <- as.data.frame(sample_data(ps_top25_genus))

# remove NA cols and cols with all identical values
V4_sample_num <- V4_sample[, sapply(V4_sample, is.numeric)]
V4_sample_num <- V4_sample_num[, !apply(is.na(V4_sample_num), 2, any)] 
V4_sample_num <- V4_sample_num[, apply(V4_sample_num, 2, function(x) length(unique(x)) > 1)] 
V4_sample_num <- V4_sample_num[, !names(V4_sample_num) %in% c("Biological_replicate", 
  "Weight_loss", "Distance_from_mine", "UXO_sum", "Main_UXO_sum", "Sum_ng_g", "Sum_ppm_ICP_MS")]

# run NMDS, extract coordinates and merge with sample data
V4_nmds <- metaMDS(V4_otu, try = 50, distance = "bray", autotransform = TRUE)
V4_nmds_df <- data.frame(MDS1 = V4_nmds$points[,1], MDS2 = V4_nmds$points[,2])
V4_nmds_df <- merge(V4_nmds_df, V4_sample, by = "row.names")
V4_nmds_df$Target <- cut(V4_nmds_df$TNT, breaks = c(-Inf, 0, Inf), labels = c("absent", "present"))
ordination_table <- V4_nmds_df

# run envfit, extract values
V4_env <- envfit(V4_nmds$points, V4_sample_num, permutations = 9999)
arrow_factor <- ordiArrowMul(V4_env)
V4_env_df <- as.data.frame(scores(V4_env, display = "vectors")) * arrow_factor
V4_env_df <- cbind(V4_env_df, Variables = rownames(V4_env_df), Pvalues = V4_env$vectors$pvals, R_squared = V4_env$vectors$r)

# select significant variables and improve labels
envfit_results_sig <- subset(V4_env_df, Pvalues < 0.001 & R_squared > 0.3)
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
extra_scaling_factor <- 2.2
colourblind <- c("#999999", "#CC79A7", "#009E73", "#D55E00", "#F0E442")
point_size <- 2.5
text_size <- 3

# plot ordination with correlating variables
nmds_envfit <-  ggplot(ordination_table, aes(x = MDS1, y = MDS2)) +
  geom_point(data = subset(ordination_table, Sample_Type == "Surface" & 
      Area %in% c("East_of_restricted", "West_of_restricted")), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(data = subset(ordination_table, Experiment != "Transect" & Area == "Restricted"), 
    aes(shape = Target), colour = "black", size = point_size + 0.75, stroke = 1.2, alpha = 0.8) +
  geom_point(aes(shape = Target, colour = Area), size = point_size, alpha = 0.7) +
  coord_fixed() +
  geom_segment(data = envfit_results_sig,
               aes(x = 0, xend = MDS1 * extra_scaling_factor * 0.95, y = 0, 
                 yend = MDS2* extra_scaling_factor * 0.95, alpha = R_squared),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_text(data = envfit_results_sig, aes(x = MDS1 * extra_scaling_factor, 
    y = MDS2 * extra_scaling_factor, label = Variables),
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
  labs(x = "nMDS1", y = "nMDS2")

ggsave(nmds_envfit, file = file.path(plot_path, "supp_nMDS_25genera_correlation.png"),
  device = "png", width = 27, height = 18, dpi = 300, unit = "cm")
