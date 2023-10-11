
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install BioPackages
BiocManager::install("BiocGenerics")
BiocManager::install("oligo")
BiocManager::install("pd.mapping50k.xba240")
BiocManager::install("hapmap100kxba")

#  Work dir change
setwd("~/Documenti/code/R/Data") 

# Load  packages
library(oligo)
library(pd.mapping50k.xba240)
library(hapmap100kxba)

if(require(pd.mapping50k.xba240) & require(hapmap100kxba))
  {
  
    celPath <- system.file("celFiles", package="hapmap100kxba")
    celFiles <- list.celfiles(celPath, full.name=TRUE)
    affySnpFeatureSet <- read.celfiles(celFiles)
  }

# test code
summary(affySnpFeatureSet)

summary(affySnpFeatureSet@assayData[["exprs"]])


