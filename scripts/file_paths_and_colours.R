#----------------- Libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DESeq2) 
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(reshape2)
library(edgeR)  
library(limma)  
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(philentropy)
library(gplots)
library(variancePartition)


#BiocManager::install("variancePartition", lib="/tgen_labs/jfryer/YOURACCOUNTHERE/R/x86_64-pc-linux-gnu-library/4.3")
#BiocManager::install("variancePartition")
#devtools::install_github("DiseaseNeuroGenomics/variancePartition")
#BiocManager::install("limma", force = TRUE) 

#----------------- Define variables
tissue <- c("Brain") # Kidney or Brain
typeOfCount <- c("ReadsPerGene.out.tab") 
#pathToRef <- c("/tgen_labs/jfryer/projects/references/mouse/ensembl_v7/")

#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

#----------------- Data
metadata <- read.delim("/tgen_labs/jfryer/kolney/chemobrain/metadata_1MPI.tsv", header = TRUE, sep = "\t")
# Update Sample_ID
#metadata$Sample_ID <- sub("MTX_LEUC", "MTXLEUC", metadata$Sample_ID)

# Exclude samples 
# metadata <- metadata[metadata$Study_Specimen_ID != "CBP_16_F_Saline",] #
