library(dada2)
library(phyloseq)
library(ggplot2)

theme_set(theme_bw())

# set up folders and files
projectPath <- "/data/projects/2019/tnt"
databasePath <- "/data/db"
plotPath <- file.path(projectPath, "analysis/plots")
savePath <- file.path(projectPath, "analysis/16s")
trainingFasta <- file.path(databasePath,
                           "silva_nr_v132_train_set.fa.gz")
speciesFasta <- file.path(databasePath,
                          "silva_species_assignment_v132.fa.gz")
                          
# load the sequence tables
seqtab1 <- readRDS(file.path(savePath, "plate1seqTable.RDS"))
seqtab2 <- readRDS(file.path(savePath, "plate2seqTable.RDS"))
seqtab3 <- readRDS(file.path(savePath, "plate3seqTable.RDS"))
seqtab4 <- readRDS(file.path(savePath, "plate4seqTable.RDS"))
stables <- list(seqtab1, seqtab2, seqtab3, seqtab4)

# merge the sequencetable from different runs into one table
all_seqTables <- mergeSequenceTables(tables = stables,
									  repeats = "error", 
									  orderBy = "abundance")                             
taxa <- assignTaxInformation(all_seqTables, trainingFasta, speciesFasta)

saveRDS(taxa, file.path(savePath, "taxonomy.RDS"))
saveRDS(all_seqtables, file.path(savePath, "merged_seqTable.RDS"))







