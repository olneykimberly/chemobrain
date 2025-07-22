#----------------- Libraries
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
library(SeuratObject)
library(Signac)
library(Seurat) 
library(stringr)
library(ggplot2)
library(harmony)
library(remaCor)
library(gridExtra)
library(grid)
library(lattice)
library(R.utils)
library(SeuratWrappers)
library(Azimuth)
library(dittoSeq)
library(dplyr)
library(RColorBrewer)
library(DESeq2) # adds matrix
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(harmony)
library(DoubletFinder)
library(reshape2)
library(ggtree)
library(BiocParallel) 
library(edgeR)  
library(limma)  
library(ggrepel) 
library(ggplot2) 
library(gplots) 
library(grDevices)  
#library(philentropy) 
library(stringr) 
library(remaCor)
library(scales)
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(data.table)
library(openxlsx)
library(gprofiler2)
library(ggpubr)
library(forcats)
library(stringr)
library(reshape2)
library(philentropy)

#BiocManager::install("variancePartition", lib="/tgen_labs/jfryer/kolney/R/x86_64-pc-linux-gnu-library/4.3")
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
metadata <- read.delim("/tgen_labs/jfryer/kolney/chemobrain/metadata.tsv", header = TRUE, sep = "\t")
# Update SampleID
metadata$SampleID <- sub("WGMRS.*", "WGMRS", metadata$SampleID)

metadata <- metadata[metadata$Study_Specimen_ID != "CBP_16_F_Saline", ]
#-----------------------
# function
