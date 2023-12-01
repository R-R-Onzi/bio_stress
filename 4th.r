devtools::install_version("dbplyr", version = "2.3.4")
BiocManager::install(version = "3.18")
BiocManager::install(gseGO)
BiocManager::install("clusterProfiler")


library(stringr)
library("DESeq2")
library("GEOquery")
library("dplyr")
library("org.Hs.eg.db")
library("conflicted")
library("clusterProfiler")

dft <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)
df <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

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


#get gene ensemble id

require("biomaRt")

conflicts_prefer(biomaRt::select)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

ens <- c("ENSG00000100601.5", "ENSG00000178826.6",
         "ENSG00000243663.1", "ENSG00000138231.8")
ensLookup <- gsub("\\.[0-9]*$", "", ens)



# ol'code

annotLookup1 <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=genes1)

annotLookup2 <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=genes2)

# new code
res_DF1 <- na.omit(res_DF1)
res_DF2 <- na.omit(res_DF2)

genes1 <- rownames(res_DF1) 
genes2 <- rownames(res_DF2) 

genes_name <- mapIds(org.Hs.eg.db, keys = genes1, column = "SYMBOL", keytype = "ENSEMBL") #map them to gene names
res_DF1$gene_name <- genes_name
res_DF1$gene_name[res_DF1$gene_name == "<NA>"] <- NA
res_DF1 <- res_DF1[!is.na(res_DF1$gene_name), ]#obtaining the genes that have a representative name

res_DF1 <-res_DF1[order(-res_DF1$log2FoldChange),]
res_DF1 <- res_DF1[res_DF1$padj<=0.05,]

genes_list1<- res_DF1$log2FoldChange
names(genes_list1) <- rownames(res_DF1)

genes_name<-mapIds(org.Hs.eg.db, keys = genes2, column = "SYMBOL", keytype = "ENSEMBL") #map them to gene names
res_DF2$gene_name <- genes_name
res_DF2$gene_name[res_DF2$gene_name == "<NA>"] <- NA
res_DF2 <- res_DF2[!is.na(res_DF2$gene_name), ]#obt

res_DF2 <-res_DF2[order(-res_DF2$log2FoldChange),]
res_DF2 <- res_DF2[res_DF2$padj<=0.05,]

genes_list2 <- res_DF2$log2FoldChange
names(genes_list2) <- rownames(res_DF2)

gse1 <- gseGO(genes_list1,  #Gene Set Enrichment Analysis
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)

gse1df <- data.frame(gse1)
write.csv(gse1df, "condition_C9_vs_C9GC.csv")


gse2 <- gseGO(genes_list2,  #Gene Set Enrichment Analysis
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)

gse2df <- data.frame(gse1)

write.csv(gse2df, "condition_C9KO_vs_C9GC.csv")


dftn <- read.csv("SG_genes.csv",header = F)
dftn <- dftn$V1
conflicts_prefer(dplyr::select)

library(hash)
## hash-2.2.6 provided by Decision Patterns
h <- hash() 
rez_list <- list()

conflicts_prefer(dplyr::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

for (ind in seq(1:length(gse1df$core_enrichment))){
  
  ensg_<-unlist(strsplit(gse1df$core_enrichment[ind],"/"))
  enrch_gene_names <- res_DF1 %>% filter(rownames(res_DF1) %in% ensg_) %>% select(gene_name)
  enrch_gene_names <- enrch_gene_names$gene_name
  if(! (1 >length(intersect(enrch_gene_names, dftn))))
    h[[gse1df$ID[ind]]] <- intersect(enrch_gene_names, dftn)
  else
    rez_list <- append(rez_list, intersect(enrch_gene_names,dftn))

}

h <- na.omit(h)


# why

conflicts_prefer(dplyr::select)

library("pheatmap")

ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

library("ggplot2")

df <- as.data.frame(colData(dds)[,c("condition","geo_accession")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
cols <- c("padj  < 0.05  & log2FoldChange > 0" = "yellow", "log2FoldChange <= 0 & padj < 0.05" = "blue", "padj > 0.05" = "red")


volcano <- ggplot(res_DF1 , aes(x =log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color =  ifelse(padj  < 0.05  & log2FoldChange > 0 , "yellow",ifelse(log2FoldChange <= 0 & padj < 0.05, "blue", "red"))), size=.001) +
  scale_color_manual(values = c("yellow", "red", "blue"), labels = c("Under expressed", "padj > 0.05", "Over expressed")) +
  xlim(-5,8) +
  ylim(0,100) +
  theme_minimal() +
  labs(
    title = "condition_C9_vs_C9GC",
    x = "Log2FoldChange",
    y = "-Log10(p-value)"
  )

print(volcano)



volcano <- ggplot(res_DF2 , aes(x =log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color =  ifelse(padj  < 0.05  & log2FoldChange > 0 , "yellow",ifelse(log2FoldChange <= 0 & padj < 0.05, "blue", "red"))), size=.001) +
  scale_color_manual(values = c("yellow", "red", "blue"), labels = c("Under expressed", "padj > 0.05", "Over expressed")) +
  xlim(-5,5) +
  ylim(0,75) +
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

filtered_genes1 <- filter(res_DF1, padj<0.05, log2FoldChange > 1 | log2FoldChange < 1)
filtered_genes2 <- filter(res_DF2, padj<0.05, log2FoldChange > 1 | log2FoldChange < 1)

resOrdered1 <- filtered_genes1[order(abs(filtered_genes1$log2FoldChange), decreasing = T),]
resOrdered2 <- filtered_genes2[order(abs(filtered_genes2$log2FoldChange), decreasing = T),]

annotLookup1 <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=row.names(resOrdered1))

annotLookup2 <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=row.names(resOrdered2))

resOrdered1$ensembl_gene_id <- row.names(resOrdered1)
resOrdered2$ensembl_gene_id <- row.names(resOrdered2)

rez1 <- resOrdered1 %>% inner_join(resOrdered1, annotLookup1, by = c("ensembl_gene_id" = "ensembl_gene_id")) 
rez2 <- resOrdered2 %>% inner_join(resOrdered2, annotLookup2, by = c("ensembl_gene_id" = "ensembl_gene_id")) 

genes_list <- resOrdered1$log2FoldChange

names(genes_list) <- rownames(genes)
genomic_idx1 <- match(resOrdered1$ensembl_gene_id, annotLookup1$ensembl_gene_id)
genomic_idx2 <- match(resOrdered2$ensembl_gene_id, annotLookup2$ensemb2_gene_id)

rpkm_ordered  <- resOrdered1[ , genomic_idx1]
rpkm_ordered  <- resOrdered2[ , genomic_idx2]

gse <- gseGO(genes_list,  #Gene Set Enrichment Analysis
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)
