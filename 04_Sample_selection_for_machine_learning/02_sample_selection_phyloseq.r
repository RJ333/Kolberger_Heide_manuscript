#!/usr/bin/env Rscript

# Select and prepare the libraries which can be used for Machine Learning
# without excessive risk of confounding variables
library(speedyseq)
library(dplyr)
library(data.table)

# select project directory, phyloseq path and the save path
projectPath <- "D:/Arbeit/projects/2019/tnt"
savePath <- file.path(projectPath, "ML/input_tables")
phyloPath <- file.path(projectPath,"phyloseq_output")
# read in phyloseq object
ps_final <- readRDS(file.path(phyloPath, "ps_final.RDS"))

# Required steps
# 1) remove unwanted samples (as prior identified)
# 2) average technical replicates,
# 3) average Mo7 (biological replicates)
# Do this for otu table and sample data and return as phyloseq object


### Sample removal
# select V4 primerset
Subset_step_1 <- subset_samples(ps_final, Primerset == "V4")

# remove libraries by station names:
stations <- c("7_GS", "9_GS", "10_GS", "ST09_60m", "ST09_100m")
Subset_step_2 <- subset_samples(Subset_step_1, !Station_ID %in% stations)

# remove libraries by Udemm numbers:
Subset_step_3 <- subset_samples(Subset_step_2, Udemm != "Udemm1157_2") 

# Get udemm numbers of MUC 2 1 slices we want to remove
core_udemms <- sample_data(subset_samples(Subset_step_3, grepl("MUC_2_1", Station_ID) & Udemm != "Udemm1326"))$Udemm
Subset_step_4 <- subset_samples(Subset_step_3, !Udemm %in% core_udemms)
# Remove all taxa columns which are now empty
Subset_step_4_filtered <- filter_taxa(Subset_step_4, function (x) {sum(x > 0) >= 1}, prune = TRUE)
# turn communities into relative abundance for averaging
Subset_relative <- transform_sample_counts(Subset_step_4_filtered, function(x) {(x / sum(x)) * 100})
otu <- otu_table(Subset_relative)
Sample_data <- as.data.frame(unclass(sample_data(Subset_relative)))
row.names(Sample_data) <- row.names(sample_data(Subset_relative))

# merge otu table with selected sample data
otu_df <- as.data.frame(as(otu, "matrix"))
Sample_data_sel <- subset(Sample_data, select = c(Udemm, Station_ID))
otu_sel <- merge(Sample_data_sel, otu_df, by = "row.names")
row.names(otu_sel) <- otu_sel$Row.names
otu_sel$Row.names <- NULL

### Libraries averaging
# calculate average per Udemm number and keep the Station_ID
setDT(otu_sel)
otu_averaged1 <- otu_sel[, lapply(.SD, mean), by = list(Udemm, Station_ID)]  # takes about 10 minutes
setDF(otu_averaged1)

# averaging ASV abundances of biological replicates of Mo7
Mo7_adjusted <- otu_averaged1 %>% filter(grepl("Mo7", Station_ID)) %>%  
  mutate(Station_ID = "Mo7", Udemm = "Udemm1350")
setDT(Mo7_adjusted)
Mo7_averaged <- Mo7_adjusted[, lapply(.SD, mean), by = list(Udemm, Station_ID)]  # about 10 minutes
setDF(Mo7_averaged)  

# removing the original, unaveraged libraries
otu_averaged2 <- otu_averaged1 %>% filter(!grepl("Mo7", Station_ID))
# adding averaged library
otu_averaged3 <- rbind(otu_averaged2, Mo7_averaged)

# add unique row names and remove sample data cols
row.names(otu_averaged3) <- otu_averaged3$Udemm
otu_averaged4 <- subset(otu_averaged3, select = c(-Udemm, -Station_ID))

### Sample data averaging, so it fits to abundance data
# prepare sample data accordingly to otu_table by removing cols which differentiate between technical replicates
drop_cols <- c("Extraction_ID", "Library", "Primerset", "Nucleic_acid", "Run", "Library_purpose", 
  "Technical_replicate", "Kit", "Kit_batch", "Extract_concentration_ng_microL")

# averaging sample data values of biological replicates of Mo7
Sample_data2 <- unique(Sample_data[, !names(Sample_data) %in% drop_cols])
Mo7_Sample_data_averaged <- Sample_data2 %>% filter(grepl("Mo7", Station_ID)) %>% 
  summarise_each(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>% 
  mutate(Station_ID = "Mo7", Udemm = "Udemm1350")

# removing the replicate entries
Sample_data3 <- Sample_data2 %>% filter(!grepl("Mo7", Station_ID))
# adding averaged entry
Sample_data4 <- rbind(Sample_data3, Mo7_Sample_data_averaged)
row.names(Sample_data4) <- Sample_data4$Udemm

# put this back into phyloseq by merging with remaining taxonomic and sequence data
Sample_data_ps <- sample_data(Sample_data4)
ps_ML0 <- phyloseq(otu_table(as.matrix(otu_averaged4), taxa_are_rows = FALSE), Sample_data_ps)
ps_ML <- merge_phyloseq(ps_ML0, refseq(ps_final), tax_table(ps_final))

# store objects
write.table(Sample_data4, file = file.path(phyloPath, "sample_data_ML.tsv"), sep = "\t")
saveRDS(ps_ML, file = file.path(phyloPath, "ps_ML.RDS"))  # selected samples and averaged replicates
