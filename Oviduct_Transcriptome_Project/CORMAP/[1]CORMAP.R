#Kenneth Lai
##setwd 
setwd("/Users/kennethlai/Desktop/HAWKINS/FAANG/RNAseq/CORMAP")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library (biomaRt)
library(limma)
library(DESeq2)
library(RColorBrewer)
library(edgeR)
#install.packages("factoextra")
library(factoextra)
library(svglite)
library(pheatmap)
library(grid)
#install.packages("gplots")
library(gplots)
#install.packages("ComplexHeatmap")
library("ComplexHeatmap")


#IMPORT DATA 
  DESeqdata=read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[global]DESEQ-results3.csv", header = T, na.strings="NA")
  phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2(edit).txt", header=T)
  phenodatareproductive= subset(phenodata, System == "Reproductive")

#Making DESeqdata more compatible 
#1a. Subsetting columns of interest (getting rid of stats info) 
  column_numbers <- c(2, 9, 10, 11, 12, 13, 14, 15, 16)
#1b. Create a new data frame with only the selected columns
  DESeqdata.filtered <- DESeqdata [, column_numbers]
#2a. Assign the first column as row names
  row.names(DESeqdata.filtered ) <- DESeqdata.filtered [, 1]

#2b. Remove the first column from the data frame
  DESeqdata.filtered  <- DESeqdata.filtered [, -1]
#DON'T NEED TO NORMALIZE? I skipped pre-PCA normalization bc DESeqdata = DESEQ2 results that were filtered by logFC >=1 and merged with already-normalized count data (ddsreproductivee)

#3a. set this DESEQ dataframe as a matrix 
  DESeqdata.filtered = as.matrix (DESeqdata.filtered)
#3b. rename all coltitles (sample names)
  colnames(DESeqdata.filtered) <- sub("^Sample_", "", colnames(DESeqdata.filtered))
  
#CREATE CORRELATION MATRIX 
cor_matrix <- cor(DESeqdata.filtered)



#COR MATRIX --> HEATMAP 


##VERSION #2 

corhmv2=Heatmap(
  cor_matrix,
  show_column_names = TRUE,
  column_title = "Correlation Matrix (DESEQ2)",
  show_row_names = TRUE,
  cluster_columns = TRUE,
  name = "Expression Correlation",
  row_title = "Samples",
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  col = colorRampPalette(c("deepskyblue3", "black","red2"))(100),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  row_names_side = "right",
  show_row_dend = TRUE,
  show_column_dend = FALSE)
  

#SAVING AS PDFs
pdf("CorrelationHeatmap(1).pdf", width = 10, height = 8)
print(corhmv2)
dev.off()




##VERSION #1 

# Set the color palette, using white for self-correlations
#my_palette <- colorRampPalette(c("deepskyblue3", "white", "red2"))(100)
#my_palette[51] <- "white"  # Set color for self-correlations to white

# Create the heatmap with adjusted color scale
#corhmv1=heatmap.2(cor_matrix,
                  #col = my_palette,
                 #main = "Correlation Map of DESeq2 Normalized Counts",
                 # xlab = "Samples",
                 # ylab = "Samples",
                 # margin = c(10, 10))


