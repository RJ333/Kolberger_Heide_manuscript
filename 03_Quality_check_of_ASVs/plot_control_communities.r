#!/usr/bin/env Rscript

# This script visualizes the ASVs present in control libraries and actual samples.
# It was furthermore used to decide which ASVs found in the controls 
# can be included in the analyses, as they likely stem from the samples rather 
# being contaminations.
# V3-V4 is not covered here as it was mostly clean and not used in the manuscript
library(speedyseq)
library(ggplot2)
library(rlang)
# select project directory, the folder containing the phyloseq objects
# and the output path for the plot
projectPath <- "D:/Arbeit/projects/2019/tnt"
phyloPath <- file.path(projectPath,"phyloseq_output")
plotPath <- file.path(projectPath, "plots")
# read in the phyloseq object
ps_full <- readRDS(file.path(phyloPath, "ps_full.RDS"))

############### NEGATIVE V4
# plot abundance of control ASVs in the rest of the samples
negative_V4_controls <- subset_samples(ps_full, 
  Library_purpose == "negative_pcr_control" & Primerset == "V4", prune = TRUE)
# filter by a threshold which will be used for machine learning anyway
negative_V4_controls <- filter_taxa(negative_V4_controls, 
  function (x) {sum(x > 0) >= 1}, prune = TRUE)      

in_neg_control <- row.names(tax_table(negative_V4_controls))
negative_V4_controls_melt <- psmelt(negative_V4_controls)
negative_V4_controls_melt$ASV <- as.factor(negative_V4_controls_melt$OTU)
negative_V4_control_ASV <- levels(droplevels(subset(negative_V4_controls_melt, 
  Abundance > 10)$ASV))
# which ASVs should be removed?
include_not <- c(
 # "ASV00001", 
  "ASV00064",
  "ASV00345",
  "ASV02051",
  "ASV02517",
  "ASV02719",
  "ASV08180",
  "ASV09367"
  )
negative_V4_control_ASV <- setdiff(negative_V4_control_ASV, include_not)
negative_V4_control_ASV2 <- row.names(tax_table(ps_full)) %in% negative_V4_control_ASV
negative_V4_control_in_samples <- prune_taxa(negative_V4_control_ASV2, ps_full)

# which Samples have the highest abundance of ASVs present in control samples?
negative_V4 <- plot_bar(negative_V4_control_in_samples, x = "Library", fill = "OTU") +
  #theme(legend.position = "none") +
  facet_wrap(~ Kit + Library_purpose, scales = "free", nrow = 2) +
  geom_text(aes(label = OTU), position = position_jitter_stack(vjust =0.5,
             jitter.height = 0.5,
             jitter.width =  0.5, offset = 1),  hjust = 0.5, size = 1.3, alpha = 0.7)
ggsave(negative_V4, filename = file.path(plotPath, "negative_V4.pdf"), device = "pdf", 
  width = 35, height = 25)

############### POSITIVE V4
# plot abundance of control ASVs in the rest of the samples
positive_V4_controls <- subset_samples(ps_full, 
  Library_purpose == "positive_pcr_control" & Primerset == "V4", prune = TRUE)
# filter by a threshold which will be used for machine learning anyway
positive_V4_controls <- filter_taxa(positive_V4_controls, 
  function (x) {sum(x > 10) >= 1}, prune = TRUE)      
positive_V4_controls_melt <- psmelt(positive_V4_controls)
positive_V4_controls_melt$ASV <- as.factor(positive_V4_controls_melt$OTU)

positive_V4_control_ASV <- levels(droplevels(subset(positive_V4_controls_melt, 
  Abundance > 10)$ASV))
include_not <- c(
  "ASV00031",
  "ASV00064",
  "ASV00087",
  "ASV00266",
  "ASV00353",
  "ASV00663",
  "ASV04593",
  "ASV06953",
  "ASV12110",
  "ASV14160",
  "ASV18038"
  )
positive_V4_control_ASV <- setdiff(positive_V4_control_ASV, include_not)
positive_V4_control_ASV2 <- row.names(tax_table(ps_full)) %in% positive_V4_control_ASV
positive_V4_control_in_samples <- prune_taxa(positive_V4_control_ASV2, ps_full)

# which Samples have the highest abundance of ASVs
positive_V4 <- plot_bar(positive_V4_control_in_samples, x = "Library", fill = "OTU") +
  #theme(legend.position = "none") +
  facet_wrap(~ Kit + Library_purpose, scales = "free", nrow = 2) +
  geom_text(aes(label = OTU), position = position_jitter_stack(vjust =0.5,
             jitter.height = 0.5,
             jitter.width =  0.5, offset = 1),  hjust = 0.5, size = 1.3, alpha = 0.7)
ggsave(positive_V4, filename = file.path(plotPath, "positive_V4_cleaned.pdf"), 
  device = "pdf", width = 35, height = 25)

############### blank V4
# plot abundance of control ASVs in the rest of the samples
blank_V4_controls <- subset_samples(ps_full, 
  Library_purpose == "blank_extraction_control" & Primerset == "V4", prune = TRUE)
# filter by a threshold which will be used for machine learning anyway
blank_V4_controls <- filter_taxa(blank_V4_controls, function (x) {sum(x > 10) >= 1}, prune = TRUE)      
blank_V4_controls_melt <- psmelt(blank_V4_controls)
blank_V4_controls_melt$ASV <- as.factor(blank_V4_controls_melt$OTU)

blank_V4_control_ASV <- levels(droplevels(subset(blank_V4_controls_melt, Abundance > 10)$ASV))
include_not <- c(
  "ASV00064",
  "ASV00345",
  "ASV02051",
  "ASV02517",
  "ASV02719",
  "ASV08180",
  "ASV09367",
  "ASV13823"
  )
blank_V4_control_ASV <- setdiff(blank_V4_control_ASV, include_not)
blank_V4_control_ASV2 <- row.names(tax_table(ps_full)) %in% blank_V4_control_ASV
blank_V4_control_in_samples <- prune_taxa(blank_V4_control_ASV2, ps_full)

# which Samples have the highest abundance of ASVs
blank_V4 <- plot_bar(blank_V4_control_in_samples, x = "Library", fill = "OTU") +
  #theme(legend.position = "none") +
  facet_wrap(~ Kit + Library_purpose, scales = "free", nrow = 2) +
  geom_text(aes(label = OTU), position = position_jitter_stack(vjust =0.5,
             jitter.height = 0.5,
             jitter.width =  0.5, offset = 1),  hjust = 0.5, size = 1.3, alpha = 0.7)
ggsave(blank_V4, filename = file.path(plotPath, "blank_V4_modified.pdf"), 
  device = "pdf", width = 35, height = 25)

################# plot control community compositions without removal of ASV
all_controls <- subset_samples(ps_full, 
  Library_purpose == "positive_pcr_control" | 
    Library_purpose == "negative_pcr_control" |
    Library_purpose == "blank_extraction_control" , prune = TRUE)
# filter by a threshold which will be used for machine learning anyway
all_controls <- filter_taxa(all_controls, function (x) {sum(x > 10) >= 1}, prune = TRUE)      
all_controls_melt <- psmelt(all_controls)
all_controls_melt$ASV <- as.factor(all_controls_melt$OTU)

all_control_ASV <- levels(droplevels(subset(all_controls_melt, Abundance > 10)$ASV))
include_not_all <- c(
  #V34
  "ASV07432",
  #V4 negative
  "ASV00345",
  "ASV02051",
  "ASV02517",
  "ASV02719",
  "ASV13897",
  # V4 positive
  "ASV00031", 
  "ASV00064",
  "ASV00087",
  "ASV00266",
  "ASV00353",
  "ASV00663",
  "ASV04593",
  "ASV06953",
  "ASV12110",
  # V4 blank
  "ASV09367",
  "ASV08180",
  "ASV13823"
  )

all_control_ASV2 <- row.names(tax_table(ps_full)) %in% all_control_ASV
all_control_in_samples <- prune_taxa(all_control_ASV2, ps_full)

# which Samples have the highest abundance of ASVs
all <- plot_bar(all_control_in_samples, x = "Library", fill = "OTU") +
  #theme(legend.position = "none") +
  facet_wrap(~ Kit + Library_purpose, scales = "free", nrow = 2) +
  geom_text(aes(label = OTU), position = position_jitter_stack(vjust =0.5,
             jitter.height = 0.5,
             jitter.width =  0.5, offset = 1),  hjust = 0.5, size = 1.3, alpha = 0.7)
ggsave(all, filename = file.path(plotPath, "all_control_ASV.pdf"), device = "pdf", 
  width = 35, height = 25)


# calculate loss of reads by removed specified ASVs
ps_full2  <- filter_taxa(ps_full, function (x) {sum(x > 10) >= 1}, prune = TRUE)      

all_ASV <- row.names(tax_table(ps_full2)) %in% include_not_all
all_in_samples <- prune_taxa(all_ASV, ps_full2)

full_library <- sample_sums(ps_full)
library_size <- sample_sums(ps_full) - sample_sums(all_in_samples)
remaining_reads <- as.data.frame(cbind(library_size, full_library))

ggplot(remaining_reads, aes(x = row.names(remaining_reads))) +
  geom_point(aes(y = full_library), colour = "black", alpha = 0.5, size = 3.5) +
  geom_point(aes(y = library_size), colour = "red", size = 1.5) +
  ggtitle("Library size after filtering and removal of selected control-present ASVs") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))
