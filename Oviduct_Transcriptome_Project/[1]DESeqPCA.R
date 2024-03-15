#Kenneth Lai
##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/PCA-DE")

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

#IMPORT DATA 
DESeqdata=read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[global]DESEQ-results7.csv", header = T, na.strings="NA")
phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2(edit).txt", header=T)
phenodatareproductive= subset(phenodata, System == "Reproductive")

#Making DESeqdata compatible w/ prcomp() for PCA analysis 
  #1a. Subetting columns of interest (getting rid of stats info) 
    column_numbers <- c(2, 9, 10, 11, 12, 13, 14, 15, 16)
  #1b. Create a new data frame with only the selected columns
    DESeqdata.filtered <- DESeqdata [, column_numbers]
  #2a. Assign the first column as row names
    row.names(DESeqdata.filtered ) <- DESeqdata.filtered [, 1]
    
  #2b. Remove the first column from the data frame
    DESeqdata.filtered  <- DESeqdata.filtered [, -1]
    #DON'T NEED TO NORMALIZE? I skipped pre-PCA normalization bc DESeqdata = DESEQ2 results that were fileterd by logFC >=1 and merged with already-normalized count data (ddsreproductivee)
    
    
    #Sample_od1_rep4 Sample_od1_rep5 Sample_od2_rep2 Sample_od2_rep6 Sample_od3_rep1 Sample_od3_rep5 Sample_ovary_rep5 Sample_ovary_rep7

    
#Perform PCA on DESeq2 results 
  #PCA-alltissues 
pdf("[6]PCA-DE_reproductive_tissues.pdf", width = 8, height = 6)
par(mar = c(8, 4, 4, 4))
par(cex = 0.8)

# Define the number of colors you want
nb.cols <- 20
#mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
pca.counts2 = prcomp(t(DESeqdata.filtered), scale=F)
#Set custom colors for each tissue
custom_colors <- c("Magnum.od1" = "#F87008", "Isthmus.od2" = "#27DAD8", "Shell_Gland.od3" = "#22E11D", "Ovary" = "#F215EE")
#Create a vector of colors based on the tissue information in phenodatareproductive$Tissue
color_vector <- custom_colors[as.character(phenodatareproductive$Tissue)]
fviz_pca_ind(pca.counts2,
             pointsize=2,
             geom.ind = ("point"), # show points only (nbut not "text")
             habillage = phenodatareproductive$Tissue, # color by groups
             palette = color_vector,
             legend.title = "Tissue",
             mean.point = FALSE
)
dev.off() 