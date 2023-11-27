library(stringr)
library("DESeq2")
library("GEOquery")
library("dplyr")

dft <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)
df <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)

meta_data <- getGEO("GSE143743")
df <- df[, -(2:6)]

row.names(df) <- df$Geneid
df <- df[, -1]

metadata_df <- pData(phenoData(meta_data[[1]]))
colnames(df)<-str_replace_all( colnames(df),"\\.", "-")
df <- df %>% select(order(colnames(df)))


df %>% select(order(colnames(df)))
metadata_df <- metadata_df[, c(1,2,10)]
metadata_df<-metadata_df[order(metadata_df$title),]

# tests
all(df$`C9KO-1_2` %in% dft$C9KO.1_2)
all(df$`C9-1_3`   %in% dft$C9.1_3)
all(df$`C9GC-1_1` %in% dft$C9GC.1_1)
all(df$`C9GC-1_3` %in% dft$C9GC.1_3)
all(df$`C9GC-1_2` %in% dft$C9GC.1_2)
all(df$`C9KO-1_3` %in% dft$C9KO.1_3)
all(df$`C9-1_2`   %in% dft$C9.1_2)
all(df$`C9-1_1`   %in% dft$C9.1_1)
all(df$`C9KO-1_1` %in% dft$C9KO.1_1)

cts <- df
coldata <- metadata_df %>%
  select(title,geo_accession)

condition <- coldata %>%
  select(title)

rownames(condition) <- coldata$title
rownames(coldata)   <- coldata$title

condition$title <- gsub("\\-.*","",condition$title)
condition$title <- gsub("\\..*","",condition$title)

coldata$title <- gsub("\\-.*","",condition$title)
coldata$title <- gsub("\\..*","",condition$title)

colnames(condition)[1] = "condition"
colnames(coldata)[1] = "condition"

condition$condition <- factor(condition$condition)
condition$condition <- relevel(condition$condition, ref="C9GC")

coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref="C9GC")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# keep only rows that have a count of at least 10 for a minimal number of samples

keep<-rowSums(counts(dds))>=10
dds<- dds[keep,]

# Differential expression analysis
dds <- DESeq(dds)

MA_plot <- plotMA(dds)


res1 <- results(dds, alpha=0.05,name="condition_C9_vs_C9GC")
res2 <- results(dds, alpha=0.05,name="condition_C9KO_vs_C9GC")

#logfolds

resultsNames(dds)

keep<-rowSums(counts(dds))>=10
dds<- dds[keep,]

resLFC1 <- lfcShrink(dds, coef="condition_C9_vs_C9GC", res = res1)

res_DF1 <- data.frame(res1)
sig_genes1 <- res_DF1 %>% filter(padj<0.05)


resLFC2 <- lfcShrink(dds, coef="condition_C9KO_vs_C9GC", res = res2)

res_DF2 <- data.frame(res2)
sig_genes2 <- res_DF2 %>% filter(padj<0.05)

# pvalues 

resOrdered1 <- res1[order(sig_genes1$log2FoldChange),]

sum(res1$padj < 0.05, na.rm=TRUE)

sig_genes1<-DESeqResults(sig_genes1)


resOrdered2 <- res2[order(sig_genes2$log2FoldChange),]

sum(res2$padj < 0.05, na.rm=TRUE)

sig_genes2<-DESeqResults(sig_genes2)

# normal
plotMA(resOrdered1, ylim=c(-2,2))

plotMA(resOrdered2, ylim=c(-2,2))


# why

library("pheatmap")

ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

library("ggplot2")

df <- as.data.frame(colData(dds)[,c("condition","geo_accession")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

volcano <- ggplot(res_DF1 , aes(x =log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color =  ifelse(padj  < 0.05  & log2FoldChange> 0 , "blue",ifelse(log2FoldChange< 0 & padj < 0.05, "red", "black")), size=.1)) +
  xlim(-5,8) +
  ylim(0,120) +
  theme_minimal() +
  labs(
    title = "condition_C9_vs_C9GC",
    x = "Log2FoldChange",
    y = "-Log10(p-value)"
  )

print(volcano)


volcano <- ggplot(res_DF2 , aes(x =log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "brown")), size = .1) +
  scale_color_manual(values = c("brown", "red")) +
  xlim(-5,5) +
  ylim(0,100) +
  theme_minimal() +
  labs(
    title = "condition_C9KO_vs_C9GC",
    x = "Log2FoldChange",
    y = "-Log10(p-value)"
  )
print(volcano)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd)

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

plotDispEsts(dds)

# map dispersion REMOVED OUTLIERS

ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                     na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)

res1 <- results(ddsCustom, alpha=0.05,name="condition_C9_vs_C9GC")
res2 <- results(ddsCustom, alpha=0.05,name="condition_C9KO_vs_C9GC")

metadata(res1)$alpha



metadata(res1)$filterThreshold

plot(metadata(res1)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res1)$lo.fit, col="red")

abline(v=metadata(res1)$filterTheta)

resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res1$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))

#logfolds

resultsNames(ddsCustom)

keep<-rowSums(counts(ddsCustom))>=10
ddsCustom<- dds[keep,]

resLFC1 <- lfcShrink(ddsCustom, coef="condition_C9_vs_C9GC", res = res1)

res_DF1 <- data.frame(res1)
sig_genes1 <- res_DF1 %>% filter(padj<0.05)


resLFC2 <- lfcShrink(ddsCustom, coef="condition_C9KO_vs_C9GC", res = res2)

res_DF2 <- data.frame(res2)
sig_genes2 <- res_DF2 %>% filter(padj<0.05)

# pvalues 

resOrdered1 <- res1[order(sig_genes1$log2FoldChange),]

sum(res1$padj < 0.05, na.rm=TRUE)

sig_genes1<-DESeqResults(sig_genes1)


resOrdered2 <- res2[order(sig_genes2$log2FoldChange),]

sum(res2$padj < 0.05, na.rm=TRUE)

sig_genes2<-DESeqResults(sig_genes2)

# normal
plotMA(resOrdered1, ylim=c(-2,2))

plotMA(resOrdered2, ylim=c(-2,2))


# why

library("pheatmap")

ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

