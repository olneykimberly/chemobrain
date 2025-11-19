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
#library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(philentropy)
library(gplots)
library(variancePartition)
library(NatParksPalettes) # colors
library(ComplexUpset)


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

### --- Function to Extract Gene Expression ---
# Helper function to extract CPM for a specific gene and merge it into the info dataframe
extract_gene_cpm <- function(gene_name, cpm_data, info_data) {
  # Subset for the specific gene
  gene_cpm <- cpm_data %>%
    dplyr::filter(gene_name == !!gene_name)
  
  # Melt and clean the data
  gene_cpm_melt <- gene_cpm %>%
    reshape2::melt() %>%
    dplyr::rename(sample = variable) %>%
    dplyr::filter(sample != "gene_name") %>% # Remove the gene_name column row
    dplyr::select(sample, value) %>%
    dplyr::rename(!!gene_name := value)
  
  # Merge the new gene column into the info dataframe
  info_data <- info_data %>%
    dplyr::left_join(gene_cpm_melt, by = "sample")
  
  return(info_data)
}

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

makePaddedDataFrame <- function(l, ...) {
  maxlen <- max(sapply(l, length))
  data.frame(lapply(l, na.pad, len = maxlen), ...)
}
#----------------- Data
#metadata_1MPI <- read.delim("/tgen_labs/jfryer/kolney/chemobrain/metadata_1MPI.tsv", header = TRUE, sep = "\t")
metadata <- read.delim("/tgen_labs/jfryer/kolney/chemobrain/metadata_1WPI.tsv", header = TRUE, sep = "\t")
metadata <- metadata %>%
  mutate(Group = gsub("MTXL", "MTXLEUC", Group)) # Why all of a sudden is it MTXL?! convert back to MTXLEUC to be consistent. 

# Update Sample_ID
# metadata$Sample_ID <- sub("MTX_LEUC", "MTXLEUC", metadata$Sample_ID)

# Exclude samples 
# "F9_SALINE", "M5_SALINE", "M3_SALINE", "M1_LEUC", "F1_SALINE", "F5_LEUC", "F4_MTX", "F1_MTX", "F6_MTXL", "M12_MTXL"
# Define the list of Sample_IDs to be removed
samples_to_remove <- c(
  "F9_SALINE", "M5_SALINE", "M3_SALINE", "M1_LEUC", "F1_SALINE",
  "F5_LEUC", "F4_MTX", "F1_MTX", "F6_MTXL", "M12_MTXL", "F4_LEUC", "M7_SALINE", "M4_MTX", "M2_MTX", "M10_MTXL", "F4_MTXL", "F8_MTXL", 
  "M11_SALINE", "M1_SALINE", "M2_SALINE", "F10_SALINE", "M11_MTXL", "F8_SALINE", "F2_SALINE", "F7_MTXL"
)

# Filter the metadata data frame to remove the specified samples 
# AND the existing condition for "CBP_16_F_Saline"
metadata <- metadata[
  !(metadata$Sample_ID %in% samples_to_remove),
]
