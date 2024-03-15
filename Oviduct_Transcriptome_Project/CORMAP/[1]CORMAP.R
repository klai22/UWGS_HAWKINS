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
library(devtools)
library(scales)
source("/Users/kennethlai/desktop/HAWKINS/FAANG/RNAseq/colorpalette/RHodor_colors/R/RHodor_palettes.R")

#IMPORT DATA 
  DESeqdata=read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[global]DESEQ-results7.csv", header = T, na.strings="NA")
  #phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2(edit).txt", header=T)
  #phenodatareproductive= subset(phenodata, System == "Reproductive")
  
  #1/23/24 : Andressa wants to rename samples w/ tissue names in figure w/ the newly abbreviated sample names in this updated file
  phenodata_updated=read.table("/Users/kennethlai/Desktop/HAWKINS/FAANG/RNAseq/CORMAP/pheno_085_final.txt", header=T)
  phenodatareproductive_updated= subset(phenodata_updated, System == "Reproductive")
  
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

#3a. Set this DESEQ dataframe as a matrix 
  DESeqdata.filtered = as.matrix (DESeqdata.filtered)
#3b. rename all coltitles (sample names)--> match to Andressa's requested sample names instead 
  colnames(DESeqdata.filtered) <- sub("^Sample_", "", colnames(DESeqdata.filtered))
  colnames(DESeqdata.filtered) <- sub("ep", "", colnames(DESeqdata.filtered))
  
#tissue names replacement? 
  #colnames(DESeqdata.filtered)[grepl("^od1", colnames(DESeqdata.filtered))] <- "Magnum(od1)"
  #colnames(DESeqdata.filtered)[grepl("^od2", colnames(DESeqdata.filtered))] <- "Isthmus(od2)"
  #colnames(DESeqdata.filtered)[grepl("^od3", colnames(DESeqdata.filtered))] <- "Shell.Gland(od3)"
  #colnames(DESeqdata.filtered)[grepl("^ovary", colnames(DESeqdata.filtered))] <- "Ovary"
  
  
#CREATE CORRELATION MATRIX 
cor_matrix <- cor(DESeqdata.filtered)



#COR MATRIX --> HEATMAP 

#DEFINING COLOR PALETTE wes anderson
#pal <- wes_palette("Zissou1", 100, type = "continuous")

#CREATING ANNOTATIONS (TISSUE - LEGEND)

  # Create a sample grouping vector
tissue_groups <- c("Magnum","Magnum", "Isthmus","Isthmus","Shell_Gland","Shell_Gland","Ovary","Ovary")

  # Convert the grouping vector to a factor
tissue_groups_factor <- factor(tissue_groups)

  # Create a color palette for groups
tissue_groups_colors <- c("Magnum" = "#F87008", "Isthmus" = "#27DAD8", "Shell_Gland" = "#22E11D", "Ovary" = "#F215EE")

  # Create the heatmap annotation based on sample groups
TissueGroupAnnot <- HeatmapAnnotation(
  df = data.frame(Tissue = tissue_groups_factor),
  col = list(Tissue = tissue_groups_colors)
)

#PLOTTING HEATMAP  
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
  #col = colorRampPalette(c("deepskyblue3", "black","red2"))(100),
  col = hodor_pal("soho") (10),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  row_names_side = "right",
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  top_annotation=TissueGroupAnnot)

#SAVING AS PDFs
pdf("CorrelationHeatmap(6).pdf", width = 10, height = 7)
print(corhmv2)
dev.off()