#set working directory
setwd("C:/Users/User/Desktop/Coding programs/R_studio/DATA_mining_R/file_project")

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("BiomaRt")

#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")

###################################################
#importing of the libraries used for the project

library(DESeq2)
library(batchtma)
library(COMBAT)
library(plyr)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(enrichR)
library(dplyr)
library(cluster)
library(dbplyr)
library(GEOquery)
library(tidyverse)
library(org.Hs.eg.db)
library(VennDiagram)
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(fuzzyjoin)
library(sf)
library(ggfortify)
library(EnhancedVolcano)
library(ggrepel)
library(conflicted)
library(clusterProfiler)
library(gseGO)
library(hash)

###################################################
# read the files containing the data and modify them if necessary


GSE128647_df <- read.table('GSE128647_htseq_raw_counts.txt', header = T, row.names = 1)

GSE234297_df <- read.table("GSE234297_gene_raw_counts.txt", header = T)

GSE235915_df <- read.table("GSE235915_genes.txt", sep = ",", row.names = 1)

GSE143743_df <- read.csv('GSE143743_counts.csv', header = T, sep = "\t", row.names = 1)

GSE195620_df <- read.table('GSE195620_U2AF1_RNA_counts.txt', header = T, row.names = 1)


################################################
#adjust the file to serve your own purpose, case for the GSE234297
#substitute Entrez gene ID with the Ensemble one 

mygenes = GSE234297_df$EntrezGeneID

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  mart = ensembl, 
  values = mygenes
)

GSE234297_df$entrezgene_id = GSE234297_df$EntrezGeneID

GSE234297_df_modified <- merge(mapping, GSE234297_df, by = "entrezgene_id")



GSE234297_df_modified = GSE234297_df_modified[ , -3]
GSE234297_df_modified = GSE234297_df_modified[ , -1]


mygenes = GSE234297_df_modified$ensembl_gene_id


row.names(GSE234297_df_modified) <- make.names(mygenes, unique = T)

GSE234297_df_modified = GSE234297_df_modified[,-1]



#remove unnecessary informations from GSE143743 and reorder 

GSE143743_df <- GSE143743_df[, -(1:5)]
colnames(GSE143743_df)<-str_replace_all( colnames(GSE143743_df),"\\.", "-")
new_order = sort(colnames(GSE143743_df))
GSE143743_df <- GSE143743_df[, new_order]
################################################


################################################
#retrieval of gene symbol for each dataset, later used for Venn Diagram and FANTOM

#GSE128647
GSE128647_gene_list <- row.names(GSE128647_df)
GSE128647_gene_list


listEnsembl()
ensembl_1 <- useEnsembl(biomart = "genes")
ens_datasets <- listDatasets(ensembl_1)

ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)


GSE128647_gene_name <- getBM( attributes = c("ensembl_gene_id", "external_gene_name" ),
                              values = GSE128647_gene_list,
                              filters = "ensembl_gene_id",
                              mart = ensembl.con)

GSE128647_gene_name



#GSE234297
GSE234297_gene_list <- row.names(GSE234297_df_modified)
GSE234297_gene_list

GSE234297_gene_name <- getBM( attributes = c("ensembl_gene_id", "external_gene_name" ),
                              values = GSE234297_gene_list,
                              filters = "ensembl_gene_id",
                              mart = ensembl.con)

GSE234297_gene_name



#GSE235915
GSE235915_gene_list <- row.names(GSE235915_df)
GSE235915_gene_list

GSE235915_gene_name <- getBM( attributes = c("ensembl_gene_id", "external_gene_name"),
                              values = GSE235915_gene_list,
                              filters = "ensembl_gene_id",
                              mart = ensembl.con)

GSE235915_gene_name



#GSE143743
GSE143743_gene_list <- row.names(GSE143743_df)
GSE143743_gene_list

GSE143743_gene_name <- getBM( attributes = c("ensembl_gene_id", "external_gene_name"),
                              values = GSE143743_gene_list,
                              filters = "ensembl_gene_id",
                              mart = ensembl.con)

GSE143743_gene_name
################################################

################################################
#check in the end the files and their features
head(GSE128647_df)
head(GSE235915_df)
head(GSE234297_df_modified)
head(GSE143743_df)

dim(GSE128647_df)
dim(GSE235915_df)
dim(GSE234297_df_modified)
dim(GSE143743_df)

class(GSE128647_df)
class(GSE235915_df)
class(GSE234297_df_modified)
class(GSE143743_df)

nrow(GSE128647_df)
nrow(GSE235915_df)
nrow(GSE234297_df_modified)
nrow(GSE143743_df)

ncol(GSE128647_df)
ncol(GSE235915_df)
ncol(GSE234297_df_modified)
ncol(GSE143743_df)

#####################################################


#####################################################

#DESeq2 normalization

#obtain the metadata for each dataframe you have
#in this case GEOquery was used due to the nature of the dataset


#GSE128647 metadata
GSE128647_metadata <- getGEO("GSE128647")
View(GSE128647_metadata)
GSE128647_metadata_df <- pData(phenoData(GSE128647_metadata[[1]]))
View(GSE128647_metadata_df)

#modify metadata 
original_values_GSE128647 <- c("GFP", "WT","FL","R151A","D247A")
replacement_values_GSE128647 <- c("Negative control", "TDP43 wild type",
                                  "2 phe mutations in RRM1","Arg mutation in RRM1 to Ala",
                                  "Asp mutation in RRM2 to Ala")


GSE128647_metadata_df$source_name_ch1 <- ifelse(GSE128647_metadata_df$source_name_ch1 %in% original_values_GSE128647,
                                                replacement_values_GSE128647[match(GSE128647_metadata_df$source_name_ch1,
                                                original_values_GSE128647)], GSE128647_metadata_df$source_name_ch1)


# choose the wanted sample information
GSE128647_metadata_df <- GSE128647_metadata_df[, c(1,2,8,19)]
View(GSE128647_metadata_df)


# Modify in order to be suitable for the study and normalization
row.names(GSE128647_metadata_df) <- GSE128647_metadata_df$title

names(GSE128647_df) <- GSE128647_metadata_df$title
View(GSE128647_df)

GSE128647_metadata_df <- GSE128647_metadata_df[, -1]





#GSE235915 metadata
GSE235915_metadata <- getGEO("GSE235915")
View(GSE235915_metadata)

GSE235915_metadata_df <- pData(phenoData(GSE235915_metadata[[1]]))
View(GSE235915_metadata_df)

# choose the wanted sample information
GSE235915_metadata_df <- GSE235915_metadata_df[, c(1,2,16,41)] 
View(GSE235915_metadata_df)

colnames(GSE235915_metadata_df)[4] <- "treatment"


row.names(GSE235915_metadata_df) <- GSE235915_metadata_df$title

names(GSE235915_df) <- GSE235915_metadata_df$title
View(GSE235915_df)

GSE235915_metadata_df <- GSE235915_metadata_df[ , -1]

#remove the known outlier from both df and metadata_df
GSE235915_df <- GSE235915_df[, -14]
GSE235915_metadata_df <- GSE235915_metadata_df[-14,] 




#GSE234297 metadata
GSE234297_metadata <- getGEO("GSE234297")
View(GSE234297_metadata)

GSE234297_metadata_df <- pData(phenoData(GSE234297_metadata[[1]]))
View(GSE234297_metadata_df)


# choose the wanted sample information
GSE234297_metadata_df <- GSE234297_metadata_df[, c(1,2,11,41,40)]  


colnames(GSE234297_metadata_df)[3] <- "Characteristics"
colnames(GSE234297_metadata_df)[4] <- "Tissue"
colnames(GSE234297_metadata_df)[5] <- "Disease_State"


row.names(GSE234297_metadata_df) <- GSE234297_metadata_df$title
names(GSE234297_df_modified) <- GSE234297_metadata_df$title
GSE234297_metadata_df <- GSE234297_metadata_df[, -1]



View(GSE234297_metadata_df)



#GSE143743 metadata
GSE143743_metadata <- getGEO("GSE143743")
View(GSE143743_metadata)

GSE143743_metadata_df <- pData(phenoData(GSE143743_metadata[[1]]))
View(GSE143743_metadata_df)

# choose the wanted sample information
GSE143743_metadata_df <- GSE143743_metadata_df[, c(1,2,10)] 
View(GSE143743_metadata_df)



row.names(GSE143743_metadata_df) <- GSE143743_metadata_df$title

names(GSE143743_df) <- GSE143743_metadata_df$title
View(GSE143743_df)

GSE143743_metadata_df <- GSE143743_metadata_df[ , -1]
colnames(GSE143743_metadata_df)[2] <- "characteristics"
###################################################
#obtain the colData variables for each set and their conditions, needed for DESeq2 application

GSE128647_df_colData <- GSE128647_metadata_df
GSE128647_df_colData$geo_accession <- as.factor(GSE128647_df_colData$geo_accession)
GSE128647_df_colData$description <- as.factor(GSE128647_df_colData$description)
GSE128647_df_colData$source_name_ch1 <- as.factor(GSE128647_df_colData$source_name_ch1)
GSE128647_df_colData$source_name_ch1 <- relevel(GSE128647_df_colData$source_name_ch1, ref = "TDP43 wild type")





GSE234297_df_colData <- GSE234297_metadata_df
GSE234297_df_colData$Disease_State <- as.factor(GSE234297_df_colData$Disease_State)
GSE234297_df_colData$Disease_State <- relevel(GSE234297_df_colData$Disease_State, ref = "Healthy control")
GSE234297_df_colData$geo_accession <- as.factor(GSE234297_df_colData$geo_accession)
GSE234297_df_colData$Characteristics <- as.factor(GSE234297_df_colData$Characteristics)
GSE234297_df_colData$Tissue <- as.factor(GSE234297_df_colData$Tissue)




GSE235915_df_colData <- GSE235915_metadata_df
GSE235915_df_colData$treatment <- as.factor(GSE235915_df_colData$treatment)
GSE235915_df_colData$geo_accession <- as.factor(GSE235915_df_colData$geo_accession)
GSE235915_df_colData$description <- as.factor(GSE235915_df_colData$description)



GSE143743_df_colData <- GSE143743_metadata_df
GSE143743_df_colData$geo_accession <- factor(GSE143743_df_colData$geo_accession)
GSE143743_df_colData$description <- factor(GSE143743_df_colData$characteristics)



# check to see if everything is how it should be (order of variables)
all(rownames(GSE128647_df_colData) %in% colnames(GSE128647_df))
all(rownames(GSE128647_df_colData) == colnames(GSE128647_df))




all(rownames(GSE234297_df_colData) %in% colnames(GSE234297_df_modified))
all(rownames(GSE234297_df_colData) == colnames(GSE234297_df_modified))




all(rownames(GSE235915_df_colData) %in% colnames(GSE235915_df))
all(rownames(GSE235915_df_colData) == colnames(GSE235915_df))




all(rownames(GSE143743_df_colData) %in% colnames(GSE143743_df))
all(rownames(GSE143743_df_colData) == colnames(GSE143743_df))

#####################################################
#DESeq2 normalization for the three datasets 


#DESeq2 for GSE128647
dds_GSE128647 <- DESeqDataSetFromMatrix(countData = GSE128647_df,
                                        colData = GSE128647_df_colData,
                                        design = ~ source_name_ch1)



#filtering
keep_val_GSE128647 <- rowSums(counts(dds_GSE128647)) >= 10
dds_GSE128647 <- dds_GSE128647[keep_val_GSE128647,]

dds_GSE128647 <- DESeq(dds_GSE128647)
res_GSE128647 <- results(dds_GSE128647,
                         contrast = c('source_name_ch1',
                                      "2 phe mutations in RRM1", "Negative control" ),
                         alpha = 0.05)

#different conditions

res_GSE128647_2 <- results(dds_GSE128647,
                         contrast = c('source_name_ch1',
                                      "Arg mutation in RRM1 to Ala", "Negative control"),
                         alpha = 0.05)


res_GSE128647_3 <- results(dds_GSE128647,
                         contrast = c('source_name_ch1',
                                      "Asp mutation in RRM2 to Ala", "Negative control"),
                         alpha = 0.05)


res_GSE128647_4 <- results(dds_GSE128647,
                           contrast = c('source_name_ch1',
                                        "TDP43 wild type", "Negative control"),
                           alpha = 0.05)



#see the results for each condition
res_GSE128647
res_GSE128647_2
res_GSE128647_3
res_GSE128647_4


summary(res_GSE128647)
summary(res_GSE128647_2)
summary(res_GSE128647_3)
summary(res_GSE128647_4)


#check the presence of NA values, if it results positive then remove them
which(is.na(res_GSE128647))
res_GSE128647 <- na.omit(res_GSE128647)
which(is.na(res_GSE128647))

which(is.na(res_GSE128647_2))
res_GSE128647_2 <- na.omit(res_GSE128647_2)
which(is.na(res_GSE128647_2))

which(is.na(res_GSE128647_3))
res_GSE128647_3 <- na.omit(res_GSE128647_3)
which(is.na(res_GSE128647_3))

which(is.na(res_GSE128647_4))
res_GSE128647_4 <- na.omit(res_GSE128647_4)
which(is.na(res_GSE128647_4))


#DESeq2 for GSE235915
dds_GSE235915 <- DESeqDataSetFromMatrix(countData = GSE235915_df,
                                        colData = GSE235915_df_colData,
                                        design = ~ treatment)


#filtering
keep_val_GSE235915 <- rowSums(counts(dds_GSE235915)) >= 10
dds_GSE235915 <- dds_GSE235915[keep_val_GSE235915,]

dds_GSE235915 <- DESeq(dds_GSE235915)
res_GSE235915 <- results(dds_GSE235915, alpha = 0.05)

res_GSE235915


#check the presence of NA values, if it results positive then remove them
which(is.na(res_GSE235915))
res_GSE235915 <- na.omit(res_GSE235915)
which(is.na(res_GSE235915))


#DESeq2 for GSE234297
dds_GSE234297 <- DESeqDataSetFromMatrix(countData = round(GSE234297_df_modified),
                                        colData = GSE234297_df_colData,
                                        design = ~ Disease_State)



#filtering
keep_val_GSE234297 <- rowSums(counts(dds_GSE234297)) >= 10
dds_GSE234297 <- dds_GSE234297[keep_val_GSE234297,]

dds_GSE234297 <- DESeq(dds_GSE234297)
res_GSE234297 <- results(dds_GSE234297, alpha = 0.05)

res_GSE235915

#check the presence of NA values, if it results positive then remove them
which(is.na(res_GSE234297))
res_GSE234297 <- na.omit(res_GSE234297)
which(is.na(res_GSE234297))


#DESeq2 for GSE143743
dds_GSE143743 <- DESeqDataSetFromMatrix(countData = GSE143743_df,
                                        colData = GSE143743_df_colData,
                                        design = ~ description)



#filtering
keep_val_GSE143743 <- rowSums(counts(dds_GSE143743)) >= 10
dds_GSE143743 <- dds_GSE143743[keep_val_GSE143743,]

#LFC shrinkage 

dds_GSE143743 <- DESeq(dds_GSE143743)
res_GSE143743 <- results(dds_GSE143743,
                         contrast = c('description',
                                      "genotype: untargetted cells containing hexanucleotide repeat expansion in intron of C9ORF72",
                                      "genotype: gene corrected (wild type C9ORF72 protein lacking the ALS mutation)"), 
                         alpha = 0.05)




res_GSE143743_2 <- results(dds_GSE143743,
                         contrast = c('description',
                                      "genotype: hexanucleotide repeat expansion in intron of C9ORF72 gene as well as a knockout of the start codon in exon 2",
                                      "genotype: gene corrected (wild type C9ORF72 protein lacking the ALS mutation)"), 
                         alpha = 0.05)


#LFC shrinkage 
resLFC1 <- lfcShrink(dds_GSE143743, res = res_GSE143743)
res_DF1 <- data.frame(res_GSE143743)
sig_genes1 <- res_GSE143743 %>% filter(padj<0.05)


resLFC2 <- lfcShrink(dds_GSE143743, res = res_GSE143743_2)
res_DF2 <- data.frame(res_GSE143743_2)
sig_genes2 <- res_GSE143743_2 %>% filter(padj<0.05)



resOrdered1 <- res_GSE143743[order(sig_genes1$log2FoldChange),]
sum(res_GSE143743$padj < 0.05, na.rm=TRUE)
sig_genes1<-DESeqResults(sig_genes1)


resOrdered2 <- res_GSE143743_2[order(sig_genes2$log2FoldChange),]
sum(res_GSE143743_2$padj < 0.05, na.rm=TRUE)
sig_genes2<-DESeqResults(sig_genes2)



#check results
res_DF1
res_DF2
#check the presence of NA values, if it results positive then remove them
which(is.na(res_DF1))
res_GSE143743 <- na.omit(res_DF1)
which(is.na(res_GSE143743))


which(is.na(res_DF2))
res_GSE143743_2 <- na.omit(res_DF2)
which(is.na(res_GSE143743_2))
###################################################


#check the common genes prior to normalisation

#this is the starting Venn diagram without normalization of the results
#first import the list of target genes that you posses and want to check
SG_related_genes <- read.delim("SG_genes.csv", sep = ",")
SG_related_genes <- SG_related_genes$ABCF1
view(SG_related_genes)

#build the struxture of the Venn Diagram
Venn_Diagram_genes_1 <- list(GSE128647 = GSE128647_gene_name$external_gene_name,
                             GSE234297 = GSE234297_gene_name$external_gene_name, 
                             GSE235915 = GSE235915_gene_name$external_gene_name,
                             GSE143743 = GSE143743_gene_name$external_gene_name,
                             SGrelated = SG_related_genes)

# create customised Venn Diagram 
ggvenn(Venn_Diagram_genes_1, show_elements = F , show_percentage = T, stroke_color="Black", 
       stroke_linetype="solid")



###################################################
#applying gene symbols to DESeq2 results


#function to obtain gene symbols and apply them to the results normalized through DESeq2
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


#function applicatoion for all scenarios

res_GSE128647$hgnc_symbol <- convertIDs( row.names(res_GSE128647), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE128647))
res_GSE128647 <- na.omit(res_GSE128647)
which(is.na(res_GSE128647))





res_GSE128647_2$hgnc_symbol <- convertIDs( row.names(res_GSE128647_2), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE128647_2))
res_GSE128647_2 <- na.omit(res_GSE128647_2)
which(is.na(res_GSE128647_2))




res_GSE128647_3$hgnc_symbol <- convertIDs( row.names(res_GSE128647_3), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE128647_3))
res_GSE128647_3 <- na.omit(res_GSE128647_3)
which(is.na(res_GSE128647_3))





res_GSE128647_4$hgnc_symbol <- convertIDs( row.names(res_GSE128647_4), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE128647_4))
res_GSE128647_4 <- na.omit(res_GSE128647_4)
which(is.na(res_GSE128647_4))


#find the common genes between the GSE128647
Venn_Diagram_genes_3 <- list(GSE128647 = res_GSE128647$hgnc_symbol,
                             GSE128647_2 = res_GSE128647_2$hgnc_symbol,
                             GSE128647_3 =res_GSE128647_3$hgnc_symbol,
                             GSE128647_4 = res_GSE128647_4$hgnc_symbol)

# create customised Venn Diagram 
ggvenn(Venn_Diagram_genes_3, show_elements = F , show_percentage = T, stroke_color="Black", 
       stroke_linetype="solid")

GSE128647_C_genes <- calculate.overlap(
  x = list(
  "GSE128647" = res_GSE128647$hgnc_symbol,
  "GSE128647_2" = res_GSE128647_2$hgnc_symbol,
  "GSE128647_3" = res_GSE128647_3$hgnc_symbol,
  "GSE128647_4" = res_GSE128647_4$hgnc_symbol
  )
)

GSE128647_C_genes <- GSE128647_C_genes[["a6"]]




res_GSE234297$hgnc_symbol <- convertIDs( row.names(res_GSE234297), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE234297))
res_GSE234297 <- na.omit(res_GSE234297)
which(is.na(res_GSE234297))





res_GSE235915$hgnc_symbol <- convertIDs( row.names(res_GSE235915), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE235915))
res_GSE235915 <- na.omit(res_GSE235915)
which(is.na(res_GSE235915))




res_GSE143743$hgnc_symbol <- convertIDs( row.names(res_GSE143743), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE143743))
res_GSE143743 <- na.omit(res_GSE143743)
which(is.na(res_GSE143743))

res_GSE143743_2$hgnc_symbol <- convertIDs( row.names(res_GSE143743_2), "ENSEMBL", "SYMBOL", org.Hs.eg.db )

#check again possible NA values
which(is.na(res_GSE143743_2))
res_GSE143743_2 <- na.omit(res_GSE143743_2)
which(is.na(res_GSE143743_2))

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

GSE143743_C_genes <- calculate.overlap(
  x = list(
    "GSE143743" = res_GSE143743$hgnc_symbol,
    "GSE143743_2" = res_GSE143743_2$hgnc_symbol
  )
)

GSE143743_C_genes <- GSE143743_C_genes[["a3"]]

Venn_Diagram_genes_4 <- list(GSE143743 = res_GSE143743$hgnc_symbol,
                             GSE143743_2 = res_GSE143743_2$hgnc_symbol)

# create customised Venn Diagram 
ggvenn(Venn_Diagram_genes_4, show_elements = F , show_percentage = T, stroke_color="Black", 
       stroke_linetype="solid")




# Venn Diagram for normalized results
Venn_Diagram_genes_2 <- list(GSE128647 = GSE128647_C_genes,
                             GSE234297 = res_GSE234297$hgnc_symbol, 
                             GSE235915 = res_GSE235915$hgnc_symbol, 
                             GSE143743 = GSE143743_C_genes,
                             SGgenes = SG_related_genes)

ggVennDiagram(Venn_Diagram_genes_2, show_elements = F , show_percentage = T,
              stroke_color="Black", label_alpha = 0,
              stroke_linetype="solid", stroke_size = 0.5,
              set_name_size = 5, label = "count") +
              ggplot2 :: scale_color_brewer(palette = "PuBu") +
              ggplot2 :: scale_fill_distiller(palette = "Blues")


common_genes <- calculate.overlap(
  x = list(
    "GSE128647" = GSE128647_C_genes,
    "GSE234297" = res_GSE234297$hgnc_symbol,
    "GSE235915" = res_GSE235915$hgnc_symbol,
    "GSE143743" = GSE143743_C_genes,
    "SGgenes" = SG_related_genes
  )
)

list_common_genes <- common_genes[["a31"]]

list_C_genes <- data.frame(list_common_genes)
write.table(list_C_genes,
            "C:/Users/User/Desktop/Coding programs/R_studio/DATA_mining_R/file_project/List_Common_genes.txt",
            quote = F, row.names = F, col.names = F)


###################################################
#MA plotting

GSE128647_MA <- plotMA(dds_GSE128647, ylim = c(-3,3))

GSE234297_MA <- plotMA(dds_GSE234297, ylim = c(-3,3))

GSE235915_MA <- plotMA(dds_GSE235915, ylim = c(-3,3))

GSE143743_MA <- plotMA(dds_GSE143743, ylim = c(-3,3))
#in the case for the GSE235915 no significant genes are presented from the MA, that is
#probably due to some errors or possibly some Batch effect coming from the dataset
#this situation will be taken in consideration again when looking again at the PCA plot

###################################################
#PCA plot

#mandatory steps before actually computing the PCA plot

vsd_GSE128647 <- vst(dds_GSE128647, blind = F)
head(assay(vsd_GSE128647), 5)

vsd_GSE234297 <- vst(dds_GSE234297, blind = F)
head(assay(vsd_GSE234297), 5)

vsd_GSE235915 <- vst(dds_GSE235915, blind = F)
head(assay(vsd_GSE235915), 5)

vsd_GSE143743 <- vst(dds_GSE143743, blind = F)
head(assay(vsd_GSE143743), 5)



#PCA_GSE128647
PCA_GSE128647 <- plotPCA(vsd_GSE128647, intgroup = c("source_name_ch1"))

PCA_GSE128647


#PCA_GSE234297
PCA_GSE234297 <- plotPCA(vsd_GSE234297, intgroup = c("Disease_State"))

PCA_GSE234297



#PCA_GSE143743
PCA_GSE143743 <- plotPCA(vsd_GSE143743, intgroup = c("description"))
PCA_GSE143743 <- PCA_GSE143743 + theme(legend.position = 'none')

PCA_GSE143743
boxplot(log10(assays(dds_GSE143743)[["cooks"]]), range=0, las=2)
plotDispEsts(dds_GSE143743)



#PCA_GSE235915
#adjust the graph and try to see if Batch effect is present
PCA_GSE235915 <- plotPCA(vsd_GSE235915, intgroup = c("treatment", "description"))
nudge <- position_nudge(y = 1)
PCA_GSE235915 + geom_text(aes(label = name), position = nudge)


#presence of Batch effect
#try to adjust the Batch effect
GSE235915_metadata_df$batch <- c("1","1","1","1","1","1","1",
                                 "2","2","2","2","2","2",
                                 "3","3","3","3","3","3","3")


#redo the normalization taking in consideration also the BATCH and factorize it
GSE235915_df_colData = GSE235915_metadata_df

GSE235915_df_colData$treatment <- as.factor(GSE235915_df_colData$treatment)
GSE235915_df_colData$geo_accession <- as.factor(GSE235915_df_colData$geo_accession)
GSE235915_df_colData$description <- as.factor(GSE235915_df_colData$description)
GSE235915_df_colData$batch <- as.factor(GSE235915_df_colData$batch)



dds_GSE235915 <- DESeqDataSetFromMatrix(countData = GSE235915_df,
                                        colData = GSE235915_df_colData,
                                        design = ~ batch+treatment)


#filtering
keep_val_GSE235915 <- rowSums(counts(dds_GSE235915)) >= 10
dds_GSE235915 <- dds_GSE235915[keep_val_GSE235915,]

dds_GSE235915 <- DESeq(dds_GSE235915)
res_GSE235915 <- results(dds_GSE235915, alpha = 0.05)

res_GSE235915


summary(res_GSE235915)

#BATCH correction
vsd_GSE235915 <- vst(dds_GSE235915, blind = F)
mat_GSE235915 <- assay(vsd_GSE235915)

mat_GSE235915 <- limma::removeBatchEffect(mat_GSE235915, vsd_GSE235915$batch)
assay(vsd_GSE235915) <- mat_GSE235915
GSE235915_df_corrected <- assay(vsd_GSE235915)


#rerun DESeq2 and the MA and PCA plot
dds_GSE235915 <- DESeqDataSetFromMatrix(countData = round(GSE235915_df_corrected),
                                        colData = GSE235915_df_colData,
                                        design = ~ treatment)


keep_val_GSE235915 <- rowSums(counts(dds_GSE235915)) >= 10
dds_GSE235915 <- dds_GSE235915[keep_val_GSE235915,]


dds_GSE235915_a <- estimateSizeFactors(dds_GSE235915)
dds_GSE235915_a <- estimateDispersionsGeneEst(dds_GSE235915_a)
dispersions(dds_GSE235915_a) <- mcols(dds_GSE235915_a)$dispGeneEst
dds_GSE235915_a <- nbinomWaldTest(dds_GSE235915_a)
res_GSE235915 <- results(dds_GSE235915_a, alpha = 0.05)
res_GSE235915

summary(res_GSE235915)

GSE235915_MA <- plotMA(dds_GSE235915, ylim = c(-3,3))

PCA_GSE235915 <- plotPCA(vsd_GSE235915, intgroup = c("treatment", "description"))
nudge <- position_nudge(y = 1)
PCA_GSE235915 + geom_text(aes(label = name), position = nudge)

#rerun DESeq2 with the contrasts that are needed for the study


#medin vs medin_ab
res_GSE235915_ab <- results(dds_GSE235915_a,
                         contrast = c('treatment',
                                      "medin_ab", "medin"),
                         alpha = 0.05)
res_GSE235915_ab
summary(res_GSE235915_ab)


#medin vs medin_nl
res_GSE235915_nl <- results(dds_GSE235915_a,
                         contrast = c('treatment',
                                      "medin_nl", "medin"),
                         alpha = 0.05)
res_GSE235915_nl
summary(res_GSE235915_nl)
#The presence of errors inside the dataset still remain, therefore it is better,
#                                       also for time reasons, to leave behind GSE235915




###################################################
#Heatmap plotting for similarities between samples

#Heatmap_GSE128647

sampledist_GSE128647 <- dist(t(assay(vsd_GSE128647)))
sampledist_matrix_GSE128647 <- as.matrix(sampledist_GSE128647)

colnames(sampledist_matrix_GSE128647) <- GSE128647_metadata_df$source_name_ch1
rownames(sampledist_matrix_GSE128647) <- GSE128647_metadata_df$source_name_ch1

colors<- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
GSE128647_Heatmap <- pheatmap(sampledist_matrix_GSE128647,clustering_distance_rows=sampledist_GSE128647,
                              clustering_distance_cols=sampledist_GSE128647,col=colors)


#Heatmap_GSE143743
ntd <- normTransform(dds_GSE143743)

select <- order(rowMeans(counts(dds_GSE143743,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df_GSE143743 <- as.data.frame(colData(dds_GSE143743)[,c("description","geo_accession")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_GSE143743)
cols <- c("padj  < 0.05  & log2FoldChange > 0" = "yellow", "log2FoldChange <= 0 & padj < 0.05" = "blue", "padj > 0.05" = "red")



#########################################################################
#Volcano plot

# Get the top 25 rows with the greatest absolute values of log2FoldChange for GSE128647
top_25_genes <- res_GSE128647_4 %>%
  arrange(desc(abs(log2FoldChange))) %>%  # Arrange by absolute log2FoldChange in descending order
  top_n(25)  # Select the top 25 genes 

pdf(file = paste0("plot_", i, ".pdf"), width = 10, height = 6)
plotCounts(dds, gene = row.names(top_25_genes)[1], intgroup = 'Condition', normalized = TRUE) # plot the most differentiated gene in each condition against all conditions
dev.off() #close pdf device
all_genes <- list()
all_genes[[i]]<-list(top25 =top_25_genes,ress=res_shrunken_df)

#do it for each condition that you take into account



#GSE128647 Volcano 




for (j in resultsNames(dds_GSE128647)){
  if (j!="Intercept"){
    
    #gene=row.names(all_genes[[j]]$res)
    
    all_genes[[j]]$ress <- all_genes[[j]]$ress %>%
      mutate(diffexp = ifelse(log2FoldChange > 0, "overexpressed", "underexpressed"))%>% filter(!is.na(padj))
    
    
    volcano_GSE128647 <-ggplot(all_genes[[j]]$ress, aes(x =log2FoldChange, y = -log10(padj),col=diffexp,label = rownames(all_genes[[j]]$ress ))) +
      
      theme_minimal() +
      labs(
        title = j,
        x = "Log2FoldChange",
        y = "-Log10(p-value)"
      )
    plot(volcano_GSE128647+geom_text(check_overlap = TRUE))
    volcano_GSE128647_2 <-ggplot(all_genes[[j]]$ress, aes(x =log2FoldChange, y = -log10(padj),col=diffexp,main=j ))+
      geom_point(aes(color = ifelse(padj  < 0.05  & diffexp == "overexpressed" , "blue",ifelse(diffexp == "underexpressed" & padj < 0.05, "red", "black"))), size = 2) +
      
      
      scale_color_manual(values = c("black", "red",'blue'),
                         labels = c("Not Significant", "Overexpressed (p < 0.05)", "Underexpressed (p < 0.05)"))
    theme_minimal() +
      labs(
        
        x = "Log2FoldChange",
        y = "-Log10(p-value)"
      )
    pdf(file = paste0("plot_volcano", j, ".pdf"), width = 10, height = 6)
    plot(volcano_GSE128647_2)
    dev.off()
  }
}



#Volcano Plot for 

volcano_GSE143743 <- ggplot(res_GSE143743 , aes(x =log2FoldChange, y = -log10(padj))) +
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

print(volcano_GSE143743)



volcano_GSE143743_2 <- ggplot(res_GSE143743_2 , aes(x =log2FoldChange, y = -log10(padj))) +
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




##########################################################

#Enrichment for the datasets

########GSE128647


for (k in resultsNames(dds_GSE128647)){
  if (k!="Intercept"){
    
    all_genes[[k]]$ress<- all_genes[[k]]$ress[order(-all_genes[[k]]$res$log2FoldChange),]
    genes<-rownames(all_genes[[k]]$res) #get gene esemble id
    
    genes_name<-mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENSEMBL") #map them to gene names
    all_genes[[k]]$ress$gene_names<- genes_name
    all_genes[[k]]$ress$gene_names[all_genes[[k]]$ress$gene_names == "<NA>"] <- NA
    all_genes[[k]]$ress<-all_genes[[k]]$ress[!is.na(all_genes[[k]]$ress$gene_names), ]#obtaining the genes that have a representative name
    
    all_genes[[k]]$ress <-all_genes[[k]]$ress[order(-all_genes[[k]]$ress$log2FoldChange),]
    all_genes[[k]]$ress <- all_genes[[k]]$ress %>%
      filter(padj <= 0.05)
    
    genes_list<- all_genes[[k]]$ress$log2FoldChange
    
    names(genes_list) <- rownames(all_genes[[k]]$ress)
    
    
    gse <- gseGO(genes_list,  #Gene Set Enrichment Analysis
                 ont = "BP",
                 keyType = "ENSEMBL",
                 OrgDb = "org.Hs.eg.db",
                 eps = 1e-300)
    all_genes[[k]]$gse<- as.data.frame(gse)
    all_genes[[k]]$gseplt<- gse
    # gseaplot(gse)
    # ridgeplot(gse)
    
    
  }
}


#Get known genes related to RNA STRESS GRANULES 
SG_genes <- read.csv("SG_genes.csv",header=F)
SG_names <- unlist(SG_genes$V1)

# for one condition: match the RNA SG with the pathways of the gene set enrichment results
common_SG <- intersect(SG_names, all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$ress$gene_names)
for (ind in seq(1:length(all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$gse$core_enrichment))){
  ensg_<-unlist(strsplit(all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$gse$core_enrichment[ind],"/"))
  enrch_gene_names<- all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$ress %>% filter(rownames(all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$ress) %in% ensg_) %>% select(gene_names)
  enrch_gene_names <- enrch_gene_names$gene_names
  
  
  
  all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$pathways[[all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$gse$Description[ind]]] <- intersect(enrch_gene_names,common_SG)
  
}
# for other condition
common_SG <- intersect(SG_names, all_genes$Condition_TDP43.wild.type_vs_Negative.control$ress$gene_names)
for (ind in seq(1:length(all_genes$Condition_TDP43.wild.type_vs_Negative.control$gse$core_enrichment))){
  
  ensg_<-unlist(strsplit(all_genes$Condition_TDP43.wild.type_vs_Negative.control$gse$core_enrichment[ind],"/"))
  enrch_gene_names<- all_genes$Condition_TDP43.wild.type_vs_Negative.control$ress %>% filter(rownames(all_genes$Condition_TDP43.wild.type_vs_Negative.control$ress) %in% ensg_) %>% select(gene_names)
  enrch_gene_names <- enrch_gene_names$gene_names
  
  
  
  all_genes$Condition_TDP43.wild.type_vs_Negative.control$pathways[[all_genes$Condition_TDP43.wild.type_vs_Negative.control$gse$Description[ind]]] <- intersect(enrch_gene_names,common_SG)
  
}


# Getting the FANTOM files of the SG that belong in the  enriched pathways
folder_path <-"C:/Users/User/Desktop/Coding programs/R_studio/DATA_mining_R/file_project"
files_in_folder <- list.files(folder_path)
genes_phantom<-c()
result_list <- list()
for (k in resultsNames(dds_GSE128647)){
  
  if (k!="Intercept" && length(all_genes[[k]]$gse!=8) ){
    selected_pathways <- names(Filter(function(x) length(x)>2, all_genes[[k]]$pathways))
    all_genes[[k]]$pathways<- all_genes[[k]]$pathways[selected_pathways] #filtering the keys
    unique_phantom_genes <- setdiff(unique(unlist( all_genes[[k]]$pathways)),genes_phantom) #SG genes in each condition uniquely
    common_SG_genes <- intersect(unique(unlist( all_genes[[k]]$pathways)), genes_phantom)
    
    
    genes_phantom <- unique(unlist(all_genes[[k]]$pathways))
    all_genes[[k]]$phantom <- unique_phantom_genes
    
    matching_files <- grep(paste("p1@(", paste(all_genes[[k]]$phantom, collapse = "|"),")(.*@.*)?\\.csv$", sep = ""), files_in_folder, value = TRUE)
    matching_files<-setdiff(matching_files,"T071118_p1@CALR3,p1@MED26.csv" )
    for (file in matching_files) { 
      data <- read.csv(paste0(folder_path, "/", file),skip=1)
      filtered_data <- data[data$Frel >= 0.98,] 
      result_list[[sub(".*@([^\\.]+)(.*@.*)?\\.csv$", "\\1",file)]]<- filtered_data
    }
  }
}




#Filtering the gse objects that match rna stress granules
all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$gseplt@result <-all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$gseplt@result %>% filter(Description %in% names(all_genes$Condition_Asp_RRM2_Ala_vs_Negative.control$pathways))

all_genes$Condition_TDP43.wild.type_vs_Negative.control$gseplt@result <-all_genes$Condition_TDP43.wild.type_vs_Negative.control$gseplt@result %>% filter(Description %in% names(all_genes$Condition_TDP43.wild.type_vs_Negative.control$pathways))



# Classifyig RNA stress granule genes found and checking which gene is common with other dataset
#classifying the the Stress granules from gene enrichment as protein or RNA


SGs <-c("RPS19", "EIF3K", "RPS6", "EIF3B", "EIF3F", "EIF3D", "PABPC1", "RPS3", "RPS24", "EIF4B",
        "EIF1", "EIF3H", "SLBP", "EIF3M", "EIF3I", "HNRNPD", "CNBP", "EIF4A1", "EIF3L", "EIF4H",
        "CSDE1", "EIF2A", "EIF3E", "HABP4", "EIF4E3", "EIF3C", "RPS11", "CLNS1A", "EIF3G", "RAN",
        "POP7", "SNRPF", "FBL", "RNASEL", "NCOA3", "MSI1", "MSI2", "FABP5", "EIF4G3", "CALR",
        "SAMD4A", "TARDBP", "RBM47", "CNOT2")

# Create a dataframe with one column "sources"
SG_df <- data.frame(sources = SGs)

RNA_vs_Protein <- read.csv("C:/Users/DELL/Desktop/RNA proj/SG_RNA.csv")
SG_df$gene_type <- ifelse(SG_df$sources %in% RNA_vs_Protein$symbol, "RNA", "Protein")

#Obtaining common SG from the two datasets
SG_from_other <- c('MKI67', 'CDC20', 'CDK1', 'CENPF', 'KIF23', 'SMC4', 'SPAG5', 'PHLDB2', 'ALPK2', 'DSP', 'TPM1', 'DTL', 'CDK2', 'BLM', 'PDLIM5', 'APOBEC3C', 'SERPINE1', 'ZFP36', 'TRIM56', 'PARP14', 'RBM47', 'YBX3', 'RBPMS', 'SOX3', 'ANXA1', 'PARK7')
intersect(SG_from_other,SGs)








#####Enrichment  GSE143743
res_GSE143743 <-res_GSE143743[order(-res_GSE143743$log2FoldChange),]
res_GSE143743 <- res_GSE143743[res_GSE143743$padj<=0.05,]

genes_list1<- res_GSE143743$log2FoldChange
names(genes_list1) <- rownames(res_GSE143743)

gse1 <- gseGO(genes_list1,  
              ont = "BP",
              keyType = "ENSEMBL",
              OrgDb = "org.Hs.eg.db",
              eps = 1e-300)

gse1df <- data.frame(gse1)

gse1df <- gse1df[abs(gse1df$NES)>=1.5,]
write.csv(gse1df, "condition_C9_vs_C9GC.csv")





res_GSE143743 <-res_GSE143743[order(-res_GSE143743$log2FoldChange),]
res_GSE143743 <- res_GSE143743[res_GSE143743$padj<=0.05,]

genes_list2 <- res_GSE143743$log2FoldChange
names(genes_list2) <- rownames(res_GSE143743)

gse2 <- gseGO(genes_list2,  
              ont = "BP",
              keyType = "ENSEMBL",
              OrgDb = "org.Hs.eg.db",
              eps = 1e-300)

gse2df <- data.frame(gse2)
gse2df <- gse2df[abs(gse2df$NES)>=1.5,]


write.csv(gse2df, "condition_C9KO_vs_C9GC.csv")

#remember to import the correct files
dftn <- read.csv("SG_genes.csv",header = F)
dftn <- dftn$V1
conflicts_prefer(dplyr::select)

library(hash)
## hash-2.2.6 provided by Decision Patterns
h <- hash() 
rez_list <- list()
rezz <- data.frame(path= "", genes= "", stringsAsFactors=FALSE)

conflicts_prefer(dplyr::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

for (ind in seq(1:length(gse1df$core_enrichment))){
  ensg_<-unlist(strsplit(gse1df$core_enrichment[ind],"/"))
  enrch_gene_names <- res_DF1 %>% filter(rownames(res_DF1) %in% ensg_) %>% select(gene_name)
  enrch_gene_names <- enrch_gene_names$gene_name
  if(! (2 >length(intersect(enrch_gene_names, dftn))))
  {
    h[[gse1df$ID[ind]]] <- intersect(enrch_gene_names, dftn)
    rezz <- rbind(rezz, list(gse1df$ID[ind], toString(intersect(enrch_gene_names, dftn))))
    rez_list <- append(rez_list, intersect(enrch_gene_names,dftn))
  }
  
}


h <- na.omit(h)

write.csv(rezz, "path_genes1.csv")

h <- hash() 
rez_list <- list()
rezz <- data.frame(path= "", genes= "", stringsAsFactors=FALSE)

for (ind in seq(1:length(gse2df$core_enrichment))){
  
  ensg_<-unlist(strsplit(gse2df$core_enrichment[ind],"/"))
  enrch_gene_names <- res_DF2 %>% filter(rownames(res_DF2) %in% ensg_) %>% select(gene_name)
  enrch_gene_names <- enrch_gene_names$gene_name
  if(! (2 >length(intersect(enrch_gene_names, dftn))))
  {
    h[[gse2df$ID[ind]]] <- intersect(enrch_gene_names, dftn)
    rezz <- rbind(rezz, list(gse2df$ID[ind], toString(intersect(enrch_gene_names, dftn))))
    rez_list <- append(rez_list, intersect(enrch_gene_names,dftn))
  }
  
}

h <- na.omit(h)

write.csv(rezz, "path_genes2.csv")

filtered_genes1 <- filter(res_GSE143743, padj<0.05, log2FoldChange > 1 | log2FoldChange < 1)
filtered_genes2 <- filter(res_GSE143743_2, padj<0.05, log2FoldChange > 1 | log2FoldChange < 1)

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

###########################################################
