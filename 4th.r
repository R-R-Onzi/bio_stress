dft <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)
df <- read.delim("4th/GSE143743_bfx1299.genecounts.csv",header=T)

meta_data <- getGEO("GSE143743")
df <- df[, -(2:6)]

row.names(df) <- df$Geneid
df <- df[, -1]

library(stringr)

metadata_df <- pData(phenoData(meta_data[[1]]))
dfd<-str_replace( colnames(df),"-", ".")
df %>% select(order(colnames(df)))



df %>% select(order(colnames(df)))
metadata_df <- metadata_df[, c(1,2,10)]
# tests
all(df$C9KO.1_2 %in% dft$C9KO.1_2)
all(df$C9.1_3   %in% dft$C9.1_3)
all(df$C9GC.1_1 %in% dft$C9GC.1_1)
all(df$C9GC.1_3 %in% dft$C9GC.1_3)
all(df$C9GC.1_2 %in% dft$C9GC.1_2)
all(df$C9KO.1_3 %in% dft$C9KO.1_3)
all(df$C9.1_2 %in% dft$C9.1_2)
all(df$C9.1_1 %in% dft$C9.1_1)
all(df$C9KO.1_1 %in% dft$C9KO.1_1)

metadata_df<-metadata_df[order(metadata_df$title), ]


all(rownames(df) %in% colnames(metadata_df))

