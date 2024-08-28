# Module Transfer Genetic Algorithm

## Description
This project is part of a masters thesis named "Transferring knowledge from common to rare diseases via a genetic algorithm".
It includes (1) R code to preprocess datasets of two diseases and to construct gene co-expression networks, (2) python code to transfer knowledge about a biological process in one disease to the other disease,
after which it is optimized to represent how the process manifests in the other disease, (3) (parts of) the data that was used, (4) python code used to analyze the results.

This algorithm can be used to identify genes that are highly related to another set of genes given as input, in the context of their co-expression and link to the disease. 

## Table of Contents
- [Requirements](#requirements)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Requirements
1. data processing (R code)
   - required libraries
     # Install BiocManager if not already installed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    # Install Bioconductor packages
    BiocManager::install(c("DESeq2", "apeglm", "WGCNA", "biomaRt"))

    # Install CRAN packages
    install.packages(c("pheatmap", "ggrepel", "gridExtra", "ggfortify", 
                   "dplyr", "limma", "flashClust", "ggplot2"))
   - instructions
     - preprocessing_alzheimer.R: uses data from Synapse (https://www.synapse.org/) for which you need to request permission.
     - preprocessing_huntington.R: uses the genes we have for alzheimer after preprocessing during the mapping of gene identifiers. As these datasets can't be shared, I uploaded a dataset called
        "background_genelist.txt" that can be used instead. Make sure to change the path from "genes_ad" to this file.
     - networks_and_clustering.R: change "soft_power" to 8 to reproduce the co-expression networks for the AD data, use power 10 to reproduce the networks for the HD data. Don't run lines 41 and 42 if you are
       running for AD.
