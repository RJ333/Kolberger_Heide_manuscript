#!/usr/bin/env Rscript

# prepare input from phyloseq object for ranger or keras classification
library(speedyseq)
library(phyloseq2ML)
library(futile.logger)
flog.threshold(INFO)

projectPath <- "/data/projects/2019/tnt"
savePath <- file.path(projectPath, "ML/input_tables")
phyloPath <- file.path(projectPath,"phyloseq_output")
ps_ML <- readRDS(file.path(phyloPath, "ps_ML.RDS"))

# set parameters
analysis_ID <- "ave_run13"
mode <- "binary_class"  # binary_class, multi_class or regression
split_ratios <- 0.75
selected_taxa <- setNames(
   c("To_Genus"), 
   c("Genus")
)

my_seed <- readr::parse_number(Sys.info()[["nodename"]])  # seed for chandler-1 is 1
copies <- 0
noise_factor <- 0.40
thresholds <- 0.08
samples <- 1
#thresholds <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1) # median of library sizes ~ 63000, 10 reads = 0.015 %

# additional sample data can be added to the community composition data:

# desired_sample_data <- c("Hg_microg_kg", "Pb206_207_ratio", "P_percent", "Ca_percent",
  # "Fe_percent", "Sc45_ppm", "V51_ppm", "Cr52_ppm", "Mn55_ppm", "Co59_ppm", "Ni60_ppm",
  # "Cu63_ppm", "Zn66_ppm", "As75_ppm", "Sr86_ppm", "Zr90_ppm", "Mo95_ppm", "Ag107_ppm", 
  # "Cd111_ppm", "Sn118_ppm", "Sb121_ppm", "Cs133_ppm", "Ba137_ppm", "W186_ppm", 
  # "Tl205_ppm", "Pb207_ppm", "Bi209_ppm", "Th232_ppm", "U238_ppm","Sum_ppm_ICP_MS", "TIC", 
  # "TN", "TC", "TS", "TOC", "Microm_001_mean", "Microm_63_mean", "Microm_125_mean", 
  # "Microm_250_mean", "Microm_500_mean", "Microm_1000_mean")

# desired_sample_data <- c("As75_ppm", "Microm_63_mean", "Microm_250_mean", 
  # "Microm_500_mean", "V51_ppm", "Zn66_ppm", "TN", "Co59_ppm", "Fe_percent")

desired_response_vars <- "TNT"

# a selection of taxa can be used as input:

# top25_genera <- c(
  # "Bacteria_Cyanobacteria_Sericytochromatia_NA_NA_NA",
  # "Bacteria_Bacteroidetes_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Maribacter",
  # "Bacteria_Proteobacteria_Gammaproteobacteria_Oceanospirillales_Halomonadaceae_Cobetia",  
  # "Bacteria_Planctomycetes_Planctomycetacia_Pirellulales_Pirellulaceae_Blastopirellula",  
  # "Bacteria_Bacteroidetes_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Maritimimonas",
  # "Bacteria_Proteobacteria_Gammaproteobacteria_Oceanospirillales_Nitrincolaceae_Motiliproteus",  
  # "Bacteria_Proteobacteria_Gammaproteobacteria_Alteromonadales_Colwelliaceae_NA",
  # "Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Thiotrichaceae_Cocleimonas",
  # "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Anaeromicrobium",
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfuromonadales_Desulfuromonadaceae_Desulfuromusa",
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_MSBL7", 
  # "Bacteria_Proteobacteria_Deltaproteobacteria_SAR324_clade(Marine_group_B)_NA_NA",
  # "Bacteria_Proteobacteria_Alphaproteobacteria_Caulobacterales_Hyphomonadaceae_NA",  
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Bdellovibrionales_Bdellovibrionaceae_OM27_clade",
  # "Bacteria_Gemmatimonadetes_Gemmatimonadetes_Gemmatimonadales_Gemmatimonadaceae_NA",
  # "Bacteria_Chloroflexi_Anaerolineae_SBR1031_A4b_NA",
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_Desulfocapsa",
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfobacterales_Desulfobulbaceae_Desulfobulbus",
  # "Bacteria_Proteobacteria_Alphaproteobacteria_Sphingomonadales_Sphingomonadaceae_Altererythrobacter",
  # "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_NA_NA",
  # "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodovibrionales_Kiloniellaceae_NA",
  # "Bacteria_Proteobacteria_Deltaproteobacteria_Desulfuromonadales_Geobacteraceae_Geopsychrobacter",    
  # "Archaea_Euryarchaeota_Methanobacteria_Methanobacteriales_Methanobacteriaceae_Methanobacterium",
  # "Bacteria_Planctomycetes_Planctomycetacia_Planctomycetales_Gimesiaceae_NA",
  # "Bacteria_TA06_NA_NA_NA_NA")
  
# best_genera <- tax_table(ps_ML)[,"To_Genus"] %in% top25_genera
# ps_top25 <- prune_taxa(best_genera, ps_ML)

# phyloseq objects as list
subset_list <- list(
 adam_V4_run13 = ps_ML
)

# subsetting and tax rank agglomeration
subset_list_tax <- create_community_table_subsets(
  subset_list = subset_list, 
  thresholds = thresholds,
  taxa_prefix = "ASV",
  num_samples = samples,
  tax_ranks = selected_taxa)
subset_list_df <- otu_table_to_df(subset_list = subset_list_tax)

# add sample data columns to the count table
# subset_list_extra <- add_sample_data(phyloseq_object = ps_ML, 
  # community_tables = subset_list_df, sample_data_names = desired_sample_data)
  
# get response variables
response_variables <- extract_response_variable(
  response_variables = desired_response_vars, phyloseq_object = ps_ML)

# discretize continuous data into two classes
responses <- categorize_response_variable(
  ML_mode = "classification", 
  response_data = response_variables, 
  my_breaks = c(-Inf, 0, Inf),
  class_labels = c("absent", "present"))

# merge the input tables with the response variables
#merged_input <- merge_input_response(subset_list_extra, responses)  # use this if you have added additional variables above
merged_input <- merge_input_response(subset_list_df, responses)

# prepare data for ranger
set.seed(my_seed)
splitted_ranger <- split_data(merged_input, split_ratios)
augmented_ranger <- augment(splitted_ranger, copies, noise_factor)

# prepare data for keras
keras_dummy <- dummify_input_tables(merged_input)
set.seed(my_seed)
splitted_keras <- split_data(keras_dummy, split_ratios)
augmented_keras <- augment(splitted_keras, copies, noise_factor)
scaled_keras <- scaling(augmented_keras)
ready_keras <- inputtables_to_keras(scaled_keras)
# new keras/tensorflow version throws error at first call --> repeat
if(!exists("ready_keras")) {
  ready_keras <- inputtables_to_keras(scaled_keras)
}

# save objects as input for machine learning
saveRDS(augmented_ranger, file = file.path(savePath, paste0("input_list_ranger_", analysis_ID, ".RDS")))
saveRDS(ready_keras, file = file.path(savePath, paste0("input_list_keras_", analysis_ID, ".RDS")))
