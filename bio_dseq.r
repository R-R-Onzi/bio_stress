if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocConductor")

install.packages("terra")
install.packages("remotes")
install.packages("raster")
install.packages("conflicted")
install.packages("gprofiler2")
remotes::install_github("vitkl/regNETcmap")
install.packages("raster")
install.packages("conflicted")
BiocManager::install("EnsDb.Hsapiens.v79")
install.packages("igraph")
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("apeglm")
BiocManager::install("ashr")
BiocManager::install("IHW")
BiocManager::install("vsn")
BiocManager::install("pheatmap")
BiocManager::install("sparseMatrixStats")
BiocManager::install("SparseArray")
BiocManager::install("DelayedMatrixStats")
BiocManager::install("clusterProfiler")

require('org.Mm.eg.db')

library("biomaRt")
library("regNETcmap")
library("igraph")
library("ggplot2")


library("apeglm")
library("ashr")
library("vsn")
library("pheatmap")
library("IHW")
library("vsn")
library("pheatmap")

library("DelayedMatrixStats")
library("clusterProfiler")
library("EnrichmentBrowser")
library("regNETcmap")

library("raster")
library("conflicted")


library("DESeq2")
library("GEOquery")
library("sparseMatrixStats")
library("SparseArray")
library("dplyr")
setwd("./") 

df <- read.delim("GSE234297_gene_raw_counts.tsv",header=T)

meta_data <- getGEO("GSE234297")

meta_data <- pData(phenoData(meta_data[[1]]))
new_df <- meta_data %>% select(title,geo_accession,characteristics_ch1.1)


new_df$title<-sub("PeripheralBlood_", "", new_df$title)

new_df$characteristics_ch1.1<-sub("disease state: ", "", new_df$characteristics_ch1.1)
colnames(new_df)[3] = "condition"

rownames(new_df)<-new_df$title 


rownames(df)<-df$EntrezGeneID 
library(org.Mm.eg.db)
library("EnsDb.Hsapiens.v79")
library("gprofiler2")
i <- 1
if (i<=46025)
{
  as.character(df$EntrezGeneID[1])
  case_when(
    class == 4 or Metro_status == 'Metro' ~ '#d62023',
    class == 4 or Metro_status == 'Non-metro' ~ '#d68182',
    class == 3 or Metro_status == 'Metro' ~ '#fc9126',
    class == 3 | Metro_status == 'Non-metro' ~ '#fcc48b',
    class == 2 | Metro_status == 'Metro' ~ '#83d921',
    class == 2 | Metro_status == 'Non-metro' ~ '#abd977',
    class == 1 
  )
}
s2 <- paste0("ENSG0000000", c(df$EntrezGeneID))
ass <- gconvert(df$EntrezGeneID,organism="hsapiens",target="ENTREZGENE",filter_na = F)$target
gs <- list(s1 = names(se), s2 = s2)
rez <- idMap(gs, org = "hsa", from = "ENSEMBL", to = "SYMBOL")
names(se) <- paste0("ENSG00000", c(df$EntrezGeneID))
a <-idMap(se, org = "hsa")

tryit <- getBM(attributes = c("ensembl_gene_id", 
                                     "hgnc_symbol"),
                      mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

cts <- df %>%
  select(-EntrezGeneID)

coldata <- new_df %>%
  select(geo_accession,condition)

condition <- new_df %>%
  select(condition)

coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref="Healthy control")

ensembl <- useMart(biomart = "www.ensembl.org", dataset="hsapiens_gene_ensembl")

banana <- getBM(attributes='hgnc_symbol', 
      filters = 'ensembl_gene_id', 
      values = charg2, 
      mart = ensembl)


dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ condition)

# keep only rows that have a count of at least 10 for a minimal number of samples

keep<-rowSums(counts(dds))>=10
dds<- dds[keep,]

# Differential expression analysis
dds <- DESeq(dds)

MA_plot <- plotMA(dds)

res <- results(dds, alpha=0.1)

#logfolds

resultsNames(dds)

keep<-rowSums(counts(dds))>=10
dds<- dds[keep,]

resLFC <- lfcShrink(dds, coef="condition_sALS_vs_Healthy.control", res = res)

res_shrunken_df <- data.frame(resLFC)
res_DF <- data.frame(res)
res_DF <- res_DF %>% filter(padj<0.1)

res_top <- res_DF %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice(1:10) 

class(res_top)

sig_genes <- res_top %>% filter(padj<0.1)

# pvalues 

resOrdered <- res[order(sig_genes$log2FoldChange),]

sum(res$padj < 0.05, na.rm=TRUE)

sig_genes<-DESeqResults(sig_genes)
# pvalues over .1

res05 <- results(dds, alpha=0.1)
res_sep <- results(dds,alpha=0.05, lfcThreshold = .7, altHypothesis="greaterAbs")
res_sep

# plots
# normal
plotMA(resOrdered, ylim=c(-2,2))
# log 
plotMA(resOrdered, ylim=c(-2,2))


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

resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult


#Plot counts

d <- plotCounts(df, gene=which.min(coldata$log2FoldChange), intgroup="condition", returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


d <- plotCounts(sig_genes, gene=which.max(sig_genes$log2FoldChange), intgroup="condition", 
                returnData=TRUE)

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

meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
# why

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("condition","geo_accession")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

volcano <-ggplot(df , aes(x =log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "brown")), size = 2) +
  scale_color_manual(values = c("brown", "red")) +
  theme_minimal() +
  labs(
    title = 6,
    x = "Log2FoldChange",
    y = "-Log10(p-value)"
  )
print(volcano)

gseKEGG()

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd)


