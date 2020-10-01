library(dada2)
library(ggplot2)
library(dplyr)

theme_set(theme_bw())

# set up folders and files
projectPath <- "/data/projects/2019/tnt"
rawPath <- file.path(projectPath, "raw_data/16s_reads")
plotPath <- file.path(projectPath, "analysis/plots")
filteredPath <- file.path(projectPath, "processed_reads/16s_reads/dada2")
savePath <- file.path(projectPath, "analysis/16s")

primer515f <- "GTGCCAGCMGCCGCGGTAA"
primer806r <- "GGACTACHVGGGTWTCTAAT"
primer341f <- "CCTACGGGNGGCWGCAG"
primer805r <- "GACTACHVGGGTATCTAATCC"

# set MiSeq and primer set 
filenames <- list.files(rawPath)
if (length(unique(sub(".*_(plate\\d)_.*", "\\1", filenames))) == 1) {
  plate <- unique(sub(".*_(plate\\d)_.*", "\\1", filenames))
}
if (length(unique(sub(".*_(V4|V34)_.*", "\\1", filenames))) == 1) {
  primerset <- unique(sub(".*_(V4|V34)_.*", "\\1", filenames))
}
if (primerset == "V34") {
  fwPrimer = primer341f
  rvPrimer = primer805r
} else if (primerset == "V4") {
  fwPrimer = primer515f
  rvPrimer = primer806r
}
stopifnot(primerset == "V34" | primerset == "V4")
stopifnot(plate == "plate1" | plate == "plate2" | plate == "plate3" | plate == "plate4" )

readPattern <- c("_R1_001.fastq", "_R2_001.fastq")
fwReads <- sort(intersect(list.files(rawPath, pattern=readPattern[1], full.names=TRUE), 
                          list.files(rawPath, pattern = primerset, full.names=TRUE)))
rvReads <- sort(intersect(list.files(rawPath, pattern=readPattern[2], full.names=TRUE), 
                          list.files(rawPath, pattern = primerset, full.names=TRUE)))
sampleNames <- sapply(strsplit(basename(fwReads), "_S"), `[`, 1)

filteredPattern <- c("_F_filt.fastq.gz", "_R_filt.fastq.gz")
fwReadsFiltered <- file.path(filteredPath,
                             paste0(sampleNames, filteredPattern[1]))
rvReadsFiltered <- file.path(filteredPath,
                             paste0(sampleNames, filteredPattern[2]))

pdf(file.path(plotPath, paste0(plate, "quality_profile_untrimmed.pdf")))
for (i in 1:length(fwReads)) {
    print(c(fwReads[i], rvReads[i]))
    print(plotQualityProfile(c(fwReads[i], rvReads[i])))
  }
dev.off()

# for finding the best parameters for filtering:

# fwReads_test <- head(fwReads, 5)
# rvReads_test <- head(rvReads, 5)
# fwReadsFiltered_test <- head(fwReadsFiltered, 5)
# rvReadsFiltered_test <- head(rvReadsFiltered, 5)

# filteredReads <- filterAndTrim(fwReads_test, fwReadsFiltered_test,
                                 # rvReads_test, rvReadsFiltered_test,
                                 # trimLeft=c(nchar(fwPrimer) + 11, 
                                            # nchar(rvPrimer)),
                                 # truncLen = c(300, 260),
                                 # maxN=0,
                                 # maxEE=c(4, 4),
                                 # truncQ=2,
                                 # rm.phix=FALSE,
                                 # compress=TRUE,
                                 # multithread=TRUE)

# checkPrimers(fwReads_test, rvReads_test,
             # fwReadsFiltered_test, rvReadsFiltered_test,
             # fwPrimer=fwPrimer,
             # rvPrimer=rvPrimer,
             # plotPath=plotPath)
  
# pdf(file.path(plotPath, paste0(plate, "quality_profile_trimmed.pdf")))
# for (i in 1:length(fwReads_test)) {
    # print(c(fwReadsFiltered_test[i], rvReadsFiltered_test[i]))
    # print(plotQualityProfile(c(fwReadsFiltered_test[i], rvReadsFiltered_test[i])))
  # }
# dev.off()

# fwError <- learnErrors(fwReadsFiltered_test, multithread=TRUE)
# rvError <- learnErrors(rvReadsFiltered_test, multithread=TRUE) 

# fwDerep <- derepFastq(fwReadsFiltered_test, verbose=FALSE)
# rvDerep <- derepFastq(rvReadsFiltered_test, verbose=FALSE)  
if(primerset == "V34") {

  filteredReads <- filterAndTrim(fwReads, fwReadsFiltered,
                                 rvReads, rvReadsFiltered,
                                 trimLeft=c(nchar(fwPrimer) + 11, 
                                            nchar(rvPrimer)),
                                 truncLen = c(300, 260),
                                 maxN=0,
                                 maxEE=c(4, 4),
                                 truncQ=2,
                                 rm.phix=FALSE,
                                 compress=TRUE,
                                 multithread=TRUE)
} else if (primerset == "V4") {
  
  filteredReads <- filterAndTrim(fwReads, fwReadsFiltered,
                               rvReads, rvReadsFiltered,
                               trimLeft=c(nchar(fwPrimer), 
                                          nchar(rvPrimer)),
                               maxN=0,
                               maxEE=c(4, 4),
                               truncQ=2,
                               rm.phix=FALSE,
                               compress=TRUE,
                               multithread=TRUE)
}
checkPrimers(fwReads, rvReads,
             fwReadsFiltered, rvReadsFiltered,
             fwPrimer=fwPrimer,
             rvPrimer=rvPrimer,
             plotPath=plotPath)
  
pdf(file.path(plotPath, paste0(plate, "quality_profile_trimmed.pdf")))
for (i in 1:length(fwReads)) {
    print(c(fwReadsFiltered[i], rvReadsFiltered[i]))
    print(plotQualityProfile(c(fwReadsFiltered[i], rvReadsFiltered[i])))
  }
dev.off()

# run Dada2
fwError <- learnErrors(fwReadsFiltered, nbases = 1e+09, multithread=TRUE)
rvError <- learnErrors(rvReadsFiltered, nbases = 1e+09, multithread=TRUE)

pdf(file.path(plotPath, paste0(plate, "error_plots.pdf")))
plotErrors(fwError, nominalQ=TRUE)
plotErrors(rvError, nominalQ=TRUE)
dev.off()

fwDerep <- derepFastq(fwReadsFiltered, verbose=FALSE)
rvDerep <- derepFastq(rvReadsFiltered, verbose=FALSE)

fwDada <- dada(fwDerep, err=fwError, multithread=TRUE)
rvDada <- dada(rvDerep, err=rvError, multithread=TRUE)

mergedReads <- mergePairs(fwDada, fwDerep, rvDada, rvDerep, maxMismatch=0)
seqTable <- makeSequenceTable(mergedReads)

if(primerset == "V34") {  
 
  # Inspect distribution of sequence lengths and keep expected sequence lengths
  # Sinclair et al, 2015 doi:10.1371/journal.pone.0116955
  table(nchar(getSequences(seqTable)))
  seqTable2 <- seqTable[,nchar(colnames(seqTable)) %in% seq(356, 430)]
 
  pdf(file.path(plotPath, paste0(plate, "sequence_length_distribution.pdf")))
  hist(nchar(getSequences(seqTable)), 
  	 main = "Full distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(270, 
                     500, 
                     by = 1)))
                     
  hist(nchar(getSequences(seqTable2)), 
  	 main = "Selected distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(356, 
                     430, 
                     by = 1)))
  dev.off()
} else if (primerset == "V4") {
  # Inspect distribution of sequence lengths and keep expected sequence lengths
  # Ziesemer et16 al, 2015 DOI: 10.1038/srep16498
  table(nchar(getSequences(seqTable)))
  seqTable2 <- seqTable[,nchar(colnames(seqTable)) %in% seq(231, 272)]
  
  pdf(file.path(plotPath, paste0(plate, "sequence_length_distribution.pdf")))
  hist(nchar(getSequences(seqTable)), 
  	 main = "Full distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(49, 
                     500, 
                     by = 1)))
  hist(nchar(getSequences(seqTable2)), 
  	 main = "Selected distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(231, 
                     272, 
                     by = 1)))
  dev.off()
}
# too many cores might cause error "sendMaster..."
# https://github.com/qiime2/q2-dada2/issues/119
seqTableNoChimera <- removeBimeraDenovo(seqTable2, method = "per-sample",
                                        multithread = 20)

# some statistics on sequence quantity and length distribution
if(primerset == "V34") {
  pdf(file.path(plotPath, paste0(plate, "sequence_length_distribution_no_bimera.pdf")))
  hist(nchar(getSequences(seqTableNoChimera)), 
  	 main = "Bimera-removed selected distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(356, 
                     430, 
                     by = 1)))
  dev.off()
} else if (primerset == "V4") {
  pdf(file.path(plotPath, paste0(plate, "sequence_length_distribution_no_bimera.pdf")))
  hist(nchar(getSequences(seqTableNoChimera)), 
  	 main = "Bimera-removed selected distribution of sequence lengths", 
  	 xlab = "Amplicon sequence length (bp)",
       breaks = (seq(231, 
                     272, 
                     by = 1)))
  dev.off()
}

chimera_abundance_loss <- 100 * (1 - sum(seqTableNoChimera)/sum(seqTable2))
chimera_unique_seq_loss <- 100 * (1 - (dim(seqTableNoChimera)/dim(seqTable2))[2])
cat("removeBimeraDenovo removed", chimera_abundance_loss, "% of counts and ",chimera_unique_seq_loss, "% of unique seqs")

trackReads(fwReads=fwReads,
           filteredReads=filteredReads,
           fwDenoised=fwDada,
           rvDenoised=rvDada,
           mergedReads=mergedReads,
           selectedLength=seqTable2,
           nonChimeraReads=seqTableNoChimera,
           sampleNames = sampleNames,
           output=file.path(plotPath, paste0(plate, "tracked_reads.pdf")))

saveRDS(seqTableNoChimera, file.path(savePath, paste0(plate, "seqTable.RDS")))

# save sessionInfo() to file
writeLines(capture.output(sessionInfo()), con = file.path(savePath, paste0(plate, "sessionInfo.txt")))
