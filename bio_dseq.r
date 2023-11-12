if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocConductor")

BiocManager::install()

BiocManager::install("DESeq2")
BiocManager::install("GEOquery")
BiocManager::install("DESeq2")

BiocManager::install("apeglm")
BiocManager::install("ashr")
BiocManager::install("IHW")
BiocManager::install("vsn")
BiocManager::install("pheatmap")



library(DESeq2)
library(GEOquery)
library(dplyr)
library(apeglm)
library(ashr)
library(vsn)
library(pheatmap)


setwd("./") 

df <- read.delim("GSE234297_gene_raw_counts.txt",header=T)

meta_data <- getGEO("GSE234297")

meta_data <- pData(phenoData(meta_data[[1]]))
new_df <- meta_data %>%
  select(title,geo_accession,characteristics_ch1.1)


new_df$title<-sub("PeripheralBlood_", "", new_df$title)

new_df$characteristics_ch1.1<-sub("disease state: ", "", new_df$characteristics_ch1.1)
colnames(new_df)[3] = "condition"

rownames(new_df)<-new_df$title 


rownames(df)<-df$EntrezGeneID 


cts <- df %>%
  select(-EntrezGeneID)

coldata <- new_df %>%
  select(geo_accession,condition)

condition <- new_df %>%
  select(condition)

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ condition)

# keep only rows that have a count of at least 10 for a minimal number of samples

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("sALS","Healthy control"))

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

#logfolds

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_Healthy.control_vs_sALS", type="apeglm")


# pvalues 

resOrdered <- res[order(res$pvalue),]

sum(res$padj < 0.1, na.rm=TRUE)
# pvalues over .1

res05 <- results(dds, alpha=0.1)


# plots
#normal
plotMA(res, ylim=c(-2,2))
# log 
plotMA(resLFC, ylim=c(-2,2))


# LFC shrinkage
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")


par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")



#IHW
# works better?
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult


#Plot counts

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# info

mcols(res)$description

# save

write.csv(as.data.frame(resOrdered), file="sALS_healthy.csv")

# transformation

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
# why
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","geo_accession")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
