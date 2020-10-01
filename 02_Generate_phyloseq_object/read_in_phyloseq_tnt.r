#!/usr/bin/env Rscript

# the initial phyloseq objects are created within
# this script, it uses some functions from the
# phyloseq2ML package

# The script requires the following dada2 output:
# the taxonomy table saved in RDS format
# the ASV table saved in RDS format
# the sample data file with matching row.names 
library(speedyseq)
library(phyloseq2ML)
library(futile.logger)
flog.threshold(INFO)

# specify paths to input files
projectPath <- "/data/projects/2019/tnt"
inputPath <- file.path(projectPath, "phyloseq_input")
outputPath <- file.path(projectPath, "phyloseq_output")

# read in files and adjust formatting
taxa <- readRDS(file.path(inputPath, "taxonomy.RDS"))
all_seqTables <- readRDS(file.path(inputPath, "merged_seqTable.RDS"))
metafile <- read.delim(file.path(inputPath, "Udemm_merged_sample_data.tsv"), 
  row.names = 1, header = TRUE, na.strings = c("", "NA"), check.names = TRUE)
row.names(metafile) <- metafile$Library
metafile2 <- sample_data(metafile)

# generate phyloseq object and adjust samples names for matching
ps_full0 <- phyloseq(otu_table(all_seqTables, taxa_are_rows = FALSE), 
               tax_table(taxa))
sample_names(ps_full0) <- stringr::str_remove(sample_names(ps_full0), "_F_filt.fastq.gz")
ps_full1 <- merge_phyloseq(ps_full0, metafile2)
ps_full2 <- standardize_phyloseq_headers(ps_full1, taxa_prefix = "ASV", use_sequences = TRUE)
ps_full <- add_unique_lineages(ps_full2)

# remove ASV present in negative PCR and blank extraction controls, remove control samples
ps_contamination_control <- subset_samples(ps_full, Library_purpose != "sample" & 
  Library_purpose != "positive_pcr_control")
ps_contamination_control <- filter_taxa(ps_contamination_control, 
  function (x) {sum(x > 0) >= 1}, prune = TRUE)

# ASV to be removed below were determined in "plot_control_communities.r" 
control_ASV <- c(
  # V3-V4
  "ASV07432",
  # V4 negative control
  "ASV00345",
  "ASV02051",
  "ASV02517",
  "ASV02719",
  "ASV13897",
  # V4 positive control
  "ASV00031", 
  "ASV00064",
  "ASV00087",
  "ASV00266",
  "ASV00353",
  "ASV00663",
  "ASV04593",
  "ASV06953",
  "ASV12110",
  # V4 blank control
  "ASV09367",
  "ASV08180",
  "ASV13823"
)
# remove control samples from data set
ps_samples <- subset_samples(ps_full, Library_purpose == "sample", prune = TRUE)
ps_samples <- filter_taxa(ps_samples,  function (x) {sum(x > 0) >= 1}, prune = TRUE)

good_ASV <- !row.names(tax_table(ps_samples)) %in% control_ASV
ps_final <- prune_taxa(good_ASV, ps_samples)

# Create taxonomy lookup and sample information for phyloseq-independent re-use
levels_tax_dictionary <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxa_vector_list <- create_taxonomy_lookup(ps_full, levels_tax_dictionary)
tax_dictionary <- tax_table(ps_full2)
full_meta_data <- as.data.frame(sample_data(ps_samples))

# store objects
saveRDS(ps_full, file = file.path(outputPath, "ps_full.RDS"))  # with control samples and control ASVs
saveRDS(ps_samples, file = file.path(outputPath, "ps_samples.RDS"))  # with control ASVs
saveRDS(ps_final, file = file.path(outputPath, "ps_final.RDS"))  # all control ASVs removed
saveRDS(taxa_vector_list, file = file.path(outputPath, "taxa_vector_list.RDS"))
write.table(tax_dictionary, file = file.path(outputPath, "tax_dictionary.tsv"), sep = "\t")
write.table(full_meta_data, file = file.path(outputPath, "sample_data_final.tsv"), sep = "\t")
