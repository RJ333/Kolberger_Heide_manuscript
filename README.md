# Kolberger_Heide_manuscript

This repository contains the code which was used to analyze microbial community compositions from the munitions dumpsite Kolberger Heide. 

The data is available under https://zenodo.org/record/4062263

The required data to generate the phyloseq object and for the ML analyses were originally stored in the folders phyloseq_input, phyloseq_output and ML_results. Due to the file upload system, they are all now in one folder, so
please adjust the respective paths in the scripts for each analysis.

## Folders

### 00 required functions

This folder contains helper functions and customized functions originating from the phyloseq2ML package. They can be sourced if needed.

### 01 Processing sequence data

This folder has the code which was used to process the 16S rRNA (gene) amplicon sequences. It leads to the generation of ASV tables and taxonomy tables.

### 02 Phyloseq object generation

The tables previously produced are combined with sample meta data to generate a phyloseq object. This object is further modified for various reasons and
optimized for machine learning.

### 03 16S Quality check

The generated ASV tables were checked with respect to contamination. It was analyzed which ASVs would occur in samples as well as in the negative (water), positive and blank extraction controls. A list of ASVs were identified which were likely actually present in the samples
and were therefore not removed from the data set, although they also appeared in control samples, but in strongly reduced abundance

### 04 Sample selection for ML

As described in the manuscript, connections between samples were investigated to identify samples which were to closely related and could therefore be identified rather by confounding variables.

### 05 Data preparation for ML

This folder contains the scripts to prepare either community composition data or just the sample data as input for RF and ANN analyses. This step was usually performed on several virtual machines (VM) in parallel.

### 06 Machine learning predictions

This folder contains the code about the actual analyses, the parameters were adjusted to the corresponding research question, the different runs were labelled e.g. as run_3. This step was usually performed on several virtual machines (VM) in parallel.
The analysis of the results is performed in the code responsible for the plotting of the data.

### 07 Machine learning robustness

This is the specific code to store the predictions over several runs and calculate an average prediction rate.

### 08 Machine learning ordination

This folder contains the code to use the proximity matrix of either supervised or unsupervised Random Forest analyses to perform PCA ordination and correlate environmental variables with these ordinations.

### 09 Machine learning variable importance

The code for Random Forest to retreive the variable importance from several runs/data splits, average and store it.

### 10 Figures and supplement

Most of the Figures and supplement of the manuscript stem directly from R, using the provided code.