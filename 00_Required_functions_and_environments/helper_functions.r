#' Creates complements of primer sequences
#'
#' Creates each possible orientation of the primer sequences (forward, reverse, complement, reverse complement).
#'
#' @param primer The primer sequence.
#'
#' @return Character vector of four elements (each possible orientation of primer sequence)
createOrients <- function(primer) {
    dna <- Biostrings::DNAString(primer)
    orients <- c(Forward=dna,
                 Complement=Biostrings::complement(dna),
                 Reverse=Biostrings::reverse(dna),
                 RevComp=Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))
}

#' Counts primer hits
#'
#' Counts hits of the primer sequences
#'
#' @param primer The primer sequence.
#' @param fn The fastq file.
#'
#' @return Number of hits for the primer sequence.
countPrimerHits <- function(primer, fn) {
    nhits <- Biostrings::vcountPattern(primer,
                                       ShortRead::sread(ShortRead::readFastq(fn)),
                                       fixed=FALSE)
    return(sum(nhits > 0))
}

#' Counts the sum of unique elements
#'
#' getN counts the sum of the unique objects in a dada or derep object.
#'
#' @param x The input object (dada or derep).
#'
#' @return Sum of unique elements.
getN <- function(x)
    sum(dada2::getUniques(x))

#' Track reads over the pipeline.
#'
#' The functions tracks the reads over the entire dada2 pipeline and creates a violin/boxplot over the amount of reads in each step (saved into the plots folder).
#'
#' @param fwReads Object of the forward unfiltered reads.
#' @param filteredReads Object of the filtered reads.
#' @param fwDenoised Object of forward denoised reads.
#' @param rvDenoised Object of reverse denoised reads.
#' @param mergedReads Object of merged reads.
#' @param selectedLength Object of reads of selected read length
#' @param nonChimeraReads Object of reads with chimera removed.
#' @param sampleNames Character vector of sample names
#' @param width Plot width in cm.
#' @param height Plot height in cm.
#' @param output Path for the barplot.
#'
#' @return None
#'
#' @export
trackReads <- function(fwReads, filteredReads, fwDenoised, rvDenoised,
                       mergedReads, selectedLength, nonChimeraReads, sampleNames,
                       width = 8.9, height = 8.9, output) {

   track <- cbind(filteredReads,
                  sapply(fwDenoised, getN),
                  sapply(rvDenoised, getN),
                  sapply(mergedReads, getN),
                  rowSums(selectedLength),
                  rowSums(nonChimeraReads))

   colnames(track) <- c("raw_input", "filtered", "fwDenoised", "rvDenoised",
                        "selected", "merged", "no-chimeras")
   rownames(track) <- sampleNames

   track <- as.data.frame(track)
   track["samplename"] <- row.names(track)
   track_melt <- reshape2::melt(track, id.vars = "samplename")

   ggplot(track_melt, aes(x = variable, y = value)) +
      geom_violin(size = 0.3, alpha = 0.9) +
      geom_boxplot(size = 0.3, width = 0.2, alpha = 0.8, outlier.alpha = 0.3, outlier.size = 0.6) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      
   ggsave(output, width = width, height = height, units = c("cm"))
}

#' Checks the sequences for primers.
#'
#' Looks for primer sequences before and after trimming and creates a barplot of the occurences.
#'
#' @param fwReads Forward read files.
#' @param rvReads Reverse read files.
#' @param fwReadsFiltered Forward filtered read files.
#' @param rvReadsFiltered Reverse filtered read files.
#' @param fwPrimer Forward primer sequence.
#' @param rvPrimer Reverse primer sequence.
#' @param width Plot width in cm.
#' @param height Plot height in cm.
#' @param plotPath Path for the output plots.
#'
#' @return None
#'
#' @export
checkPrimers <- function(fwReads, rvReads,
                         fwReadsFiltered, rvReadsFiltered,
                         fwPrimer, rvPrimer,
                         width = 8.9, height = 8.9,
                         plotPath) {
    fwOrients <- createOrients(fwPrimer)
    rvOrients <- createOrients(rvPrimer)

    primerBefore <- as.data.frame(
        rbind(fwForwardReads=sapply(fwOrients, countPrimerHits, fn=fwReads[[1]]),
              fwReverseReads=sapply(fwOrients, countPrimerHits, fn=rvReads[[1]]),
              rvForwardReads=sapply(rvOrients, countPrimerHits, fn=fwReads[[1]]),
              rvReverseReads=sapply(rvOrients, countPrimerHits, fn=rvReads[[1]])))
    primerBefore["id"] <- row.names(primerBefore)
    primerBefore["stage"] <- "before"

    primerAfter <- as.data.frame(
        rbind(fwForwardReads=sapply(fwOrients, countPrimerHits, fn=fwReadsFiltered[[1]]),
              fwReverseReads=sapply(fwOrients, countPrimerHits, fn=rvReadsFiltered[[1]]),
              rvForwardReads=sapply(rvOrients, countPrimerHits, fn=fwReadsFiltered[[1]]),
              rvReverseReads=sapply(rvOrients, countPrimerHits, fn=rvReadsFiltered[[1]])))
    primerAfter["id"] <- row.names(primerAfter)
    primerAfter["stage"] <- "after"

    df <- rbind(primerBefore, primerAfter)
    df <- reshape2::melt(df)

    ggplot2::ggplot(df, ggplot2::aes_string(x="id", y="value", fill="stage")) +
        ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(), colour="black") +
        ggplot2::facet_grid(. ~ variable) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.5))

    ggplot2::ggsave(file.path(plotPath, paste0(plate, "primer_presence.pdf")),
                    width = width, height = height, units = c("cm"))
}

#' Adding taxonomic information to sequence table.
#'
#' This function adds the taxonomic information from the trainings and species fasta file.
#'
#' @param seqTableNoChimera Sequence table.
#' @param trainingFasta Path to training fasta file.
#' @param speciesFasta Path to species fasta file.
#' @param ... additional parameters for dada2::assignTaxonomy or dada2::addSpecies.
#'
#' @return Sequence table with assigned taxonomy.
#'
#' @export
assignTaxInformation <- function(seqTableNoChimera, trainingFasta, speciesFasta = NULL, ...) {
    taxa <- dada2::assignTaxonomy(seqTableNoChimera, paste0(trainingFasta),
                           multithread=TRUE, ...)
    if (!is.na(speciesFasta))
        dada2::addSpecies(taxa, paste0(speciesFasta), ...)
}

#' Calculate univariate outliers as in a boxplot
#'
#' This function replaces outlier values in a data.frame column with NA
#'
#' @param x a numerical data.frame selection
#'
#' @return
#'
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#' Calculate confidence intervals
#'
#' This function calculates confidence intervals for different probabilities
#'
#' @param x a numerical vector with values
#' @param y the quantile probability as decimal e.g. 0.25 
#'
#' @return the value of the specified probability
#'
#' @export
#' @examples 
#' confidence_interval(1:1000, 0.95)
#' confidence_interval(-1000:1000, 0.10)
confidence_interval <- function(x, y) {
  quantile(x, probs = y, na.rm = TRUE)
}

#' Permute a data frame's column
#'
#' This function allows to permute selected columns of the specified data frame 
#' and returns a list of data frames with the permuted columns modified from 
#' [here](https://community.rstudio.com/t/modelr-permutate-dataframe-output/27762/7)
#'
#' @param data the dataset with the column to permute
#' @param n the number of permutations
#' @param perm_cols character vector with the columns names to be permuted
#'
#' @return A list of data frames with the permuted columns
#'
#' @export
#' @examples 
#' shuffle(data = iris, n = 10, perm_cols = "Species")
shuffle <- function(data, n, perm_cols){
  
  cols_ids <- match(perm_cols, colnames(data))
  ids <- seq_len(nrow(data))
  n_ids <- purrr::rerun(n, sample(ids))
  
  purrr::map(n_ids, function(x){
    data[ids, cols_ids] <- data[x, cols_ids]
    data
  })
}

#' position jitter for labels on stacked bars
#'
#' This function allows to combine jitter and stack in ggplot2 objects
#'
#' @param vjust vertical adjustment of labels, default = 1
#' @param reverse the orientation, default = FALSE
#' @param jitter.width the width of the jittering, default = 1
#' @param jitter.height the hight of the jittering, default = 1
#' @param jitter.seed provide a reproducible seed? Default = NULL
#' @param offset provide offset? Default = NULL
#'
#' @return a ggplot2 ggproto object
#'
#' @export
#' taken from here: https://stackoverflow.com/questions/50059193/jitter-text-labels-with-position-stack
position_jitter_stack <- function(vjust = 1, reverse = FALSE, 
                                  jitter.width = 1, jitter.height = 1,
                                  jitter.seed = NULL, offset = NULL) {
  ggproto(NULL, PositionJitterStack, vjust = vjust, reverse = reverse, 
          jitter.width = jitter.width, jitter.height = jitter.height,
          jitter.seed = jitter.seed, offset = offset)
}

PositionJitterStack <- ggproto("PositionJitterStack", PositionStack,
  type = NULL,
  vjust = 1,
  fill = FALSE,
  reverse = FALSE,
  jitter.height = 1,
  jitter.width = 1,
  jitter.seed = NULL,
  offset = 1,

  setup_params = function(self, data) {
    list(
      var = self$var %||% ggplot2:::stack_var(data),
      fill = self$fill,
      vjust = self$vjust,
      reverse = self$reverse,
      jitter.height = self$jitter.height,
      jitter.width = self$jitter.width,
      jitter.seed = self$jitter.seed,
      offset = self$offset
    )
  },

  setup_data = function(self, data, params) {
    data <- PositionStack$setup_data(data, params)
    if (!is.null(params$offset)) {
      data$to_jitter <- sapply(seq(nrow(data)), function(i) {
        any(abs(data$y[-i] - data$y[i]) <= params$offset)
      })
    } else {
      data$to_jitter <- TRUE
      }
    data
  },

  compute_panel = function(data, params, scales) {
    data <- PositionStack$compute_panel(data, params, scales)

    jitter_df <- data.frame(width = params$jitter.width,
                            height = params$jitter.height)

    if (!is.null(params$jitter.seed)) jitter_df$seed = params$jitter.seed
    jitter_positions <- PositionJitter$compute_layer(
      data[data$to_jitter, c("x", "y")],
      jitter_df
    )

    data$x[data$to_jitter] <- jitter_positions$x
    data$y[data$to_jitter] <- jitter_positions$y

    data
  }
)

#' Extract proximity matrix from ranger random forest
#'
#' This function extracts the proximity matrix of a random forest generated
#' in ranger. The proximity matrix was calculated only using out of bag (oob)
#' samples. It can be used for unsupervised learning. Taken from 
#' [here](https://github.com/imbs-hl/ranger/issues/514)
#'
#' @param fit the trained ranger object
#' @param olddata the data fit was trained upon
#'
#' @return The proximity matrix
#'
#' @export
extract_proximity_oob = function(fit, olddata) {
  pred = predict(fit, olddata, type = "terminalNodes")$predictions
  prox = matrix(NA, nrow(pred), nrow(pred))
  ntree = ncol(pred)
  n = nrow(prox)
  
  if (is.null(fit$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }
  
  # Get inbag counts
  inbag = simplify2array(fit$inbag.counts)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Use only trees where both obs are OOB
      tree_idx = inbag[i, ] == 0 & inbag[j, ] == 0
      prox[i, j] = sum(pred[i, tree_idx] == pred[j, tree_idx]) / sum(tree_idx)
    }
  }
  
  prox
}
