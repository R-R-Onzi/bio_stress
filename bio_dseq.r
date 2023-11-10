if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("GEOquery")

library(DESeq2)
library(GEOquery)
library(dplyr)

setwd("~/Documenti/code/python/bio_stress") 

df <- read.delim("GSE234297_gene_raw_counts.tsv",header=T)

meta_data <- getGEO("GSE234297")

meta_data <- pData(phenoData(meta_data[[1]]))
new_df <- meta_data %>%
  select(title,geo_accession,characteristics_ch1.1)

rownames(new_df)<-meta_data$title 
rownames(new_df)<-sub("PeripheralBlood_", "", rownames(new_df))

new_df <- new_df %>%
  select(geo_accession,characteristics_ch1.1)

new_df <- new_df %>%
  new_df$characteristics_ch1.1

rownames(df)<-df$EntrezGeneID 


df <- df %>%
  select(-EntrezGeneID)
