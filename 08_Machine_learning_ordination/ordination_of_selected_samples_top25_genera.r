library(speedyseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ranger)

projectPath <- "/data/projects/2019/tnt"
#projectPath <- "D:/Arbeit"
savePath <- file.path(projectPath, "ML/input_tables")
phyloPath <- file.path(projectPath,"phyloseq_output")
resultPath <- file.path(projectPath, "ML/results")
ML <- "ranger"
nodename <- Sys.info()[["nodename"]]
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))

num.trees <- 10000
mtry_factor <- 1

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

# reduce data set
V4_sample_num <- V4_sample[, sapply(V4_sample, is.numeric)]
V4_sample_num <- V4_sample_num[, !apply(is.na(V4_sample_num), 2, any)] # remove NA cols
V4_sample_num <- V4_sample_num[, apply(V4_sample_num, 2, function(x) length(unique(x)) > 1)] # remove cols with all identical values
V4_sample_num <- V4_sample_num[, !names(V4_sample_num) %in% c("Biological_replicate", "Weight_loss", "Distance_from_mine", "UXO_sum", "Main_UXO_sum", "Sum_ng_g", "Sum_ppm_ICP_MS")]
cols <- c("Area", "Collection", "Date", "Experiment")
V4_sample_combined <- cbind(V4_sample[,cols], V4_sample_num)
# PCA
unsupervised_table <- V4_otu
all_vars <- ncol(unsupervised_table)
for_mtry <- ifelse((sqrt(all_vars) * mtry_factor) < all_vars,
  sqrt(all_vars) * mtry_factor, all_vars)

synth_data <- as.data.frame(lapply(as.data.frame(unsupervised_table), function(x) {
  sample(x, length(x), replace = TRUE)
}))
combined_data <- rbind(data.frame(y = 0, unsupervised_table), 
             data.frame(y = 1, synth_data))
combined_data$y <- factor(combined_data$y)

# Run unsupervised Random Forest
forest <- ranger(y ~ ., combined_data, keep.inbag = TRUE, num.trees = num.trees, mtry = for_mtry)
prox <- extract_proximity_oob(forest, combined_data)[1:nrow(unsupervised_table), 1:nrow(unsupervised_table)]
row.names(prox) <- row.names(unsupervised_table)

# Generate dist matrix and % of variation for x and y axis
distance_matrix <- dist(1 - prox)
PCA_object <- cmdscale(distance_matrix, eig = TRUE, x.ret = TRUE)

#V4_env <- envfit(PCA_object, V4_sample_num, permutations = 9999)
V4_env <- envfit(PCA_object, V4_sample_combined, permutations = 9999)
arrow_factor <- ordiArrowMul(V4_env)
V4_env_df <- as.data.frame(scores(V4_env, display = "vectors")) * arrow_factor
V4_env_df <- cbind(V4_env_df, Variables = rownames(V4_env_df), Pvalues = V4_env$vectors$pvals, R_squared = V4_env$vectors$r)

PCA_variation <- round(PCA_object$eig/sum(PCA_object$eig) * 100, 1)
PCA_values <- scores(PCA_object)

# Combine to plot ready data frame
PCA_data <- data.frame(
  Sample = rownames(PCA_values),
  X = PCA_values[, 1],
  Y = PCA_values[, 2])

# merge with sample data
RF_PCA_meta <- merge(PCA_data, V4_sample, 
  by.x = "Sample", by.y = "row.names", all.x = TRUE)
RF_PCA_meta$Target <- cut(RF_PCA_meta$TNT, breaks = c(-Inf, 0, Inf), labels = c("absent", "present"))
#RF_PCA_meta$Target <- cut(RF_PCA_meta$Distance_from_mine, breaks = c(-Inf, 2, Inf), labels = c("near", "remote"))
RF_PCA_meta$Nodename <- nodename

write.table(RF_PCA_meta, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", nodename, "_URF_PCA_genus25_008.tsv")))
write.table(V4_env_df, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", nodename, "_URF_PCA_genus25_008_envdata.tsv")))
write.table(PCA_variation, sep = "\t", file = file.path(resultPath, 
  paste0(ML, "_", nodename, "_URF_PCA_genus25_008_variation.tsv")))


### plotting V4
extra_scaling_factor <- 2
ggplot(RF_PCA_meta, aes(x = X, y = Y)) +
  geom_point(aes(colour = Area, shape = Target), size = 2.5, alpha = 0.5) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = V4_env_df2,
               aes(x = 0, xend = Dim1 * extra_scaling_factor, y = 0, yend = Dim2* extra_scaling_factor),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = V4_env_df2, aes(x = Dim1 * extra_scaling_factor, y = Dim2 * extra_scaling_factor, label = Variables),
            size = 2, alpha = 0.5) +
  geom_text(aes(label = Station_ID), colour = "black", size = 1, alpha = 0.5) +
  labs(x = paste("PC1 - ", PCA_variation[1], "%", sep = ""),
    y = paste("PC2 - ", PCA_variation[2], "%", sep = ""))

ggplot(RF_PCA_meta, aes(x = X, y = Y)) +
  geom_point(aes(colour = Cruise_ID, fill = Target), stroke = 1.5, pch = 21, size = 3.5, alpha = 0.8) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_text(aes(label = Station_ID), colour = "black", size = 2, alpha = 0.5) +
  labs(x = paste("PC1 - ", PCA_variation[1], "%", sep = ""),
    y = paste("PC2 - ", PCA_variation[2], "%", sep = ""))

# run NMDS, extract coordinates and merge with sample data
V4_nmds <- metaMDS(V4_otu, try = 30, distance = "bray", autotransform = TRUE)
V4_nmds_df <- data.frame(MDS1 = V4_nmds$points[,1], MDS2 = V4_nmds$points[,2])
V4_nmds_df <- merge(V4_nmds_df, V4_sample, by = "row.names")
V4_nmds_df$Target <- cut(V4_nmds_df$TNT, breaks = c(-Inf, 0, Inf), labels = c("absent", "present"))

# run envfit, extract values and select significant variables
V4_env <- envfit(V4_nmds$points, V4_sample_num, permutations = 9999)
arrow_factor <- ordiArrowMul(V4_env)
V4_env_df <- as.data.frame(scores(V4_env, display = "vectors")) * arrow_factor
V4_env_df <- cbind(V4_env_df, Variables = rownames(V4_env_df), Pvalues = V4_env$vectors$pvals, R_squared = V4_env$vectors$r)
V4_env_df2 <- subset(V4_env_df, Pvalues < 0.001)

### plotting V4

# plot ordination with envfit arrows for numerical variables
ggplot(V4_nmds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = Target, shape = Area), stroke = 1.5, size = 3.5, alpha = 0.8) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = V4_env_df2,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(aes(label = Station_ID), size = 1.5, colour = "black") +             
  geom_text(data = V4_env_df2, aes(x = MDS1, y = MDS2, label = Variables),
            size = 3, alpha = 0.5)
