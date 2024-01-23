#Kenneth Lai
##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/DEFigures")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library(limma)
library(DESeq2)
library(RColorBrewer)
library(edgeR)
#install.packages("factoextra")
library(factoextra)
library(svglite)
#install.packages("tidyverse")
library(tidyverse)
library("pheatmap")

#IMPORT DATA - I believe that the #s for each tissue are "normalized counts" by DESEQ's internal normalization method, need to ask Andressa???
  DESEQdata= read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/DESEQ-results.csv", header=T)

  ##CLEANING DATA (isolating only genes w/ padj = or < 0.01 , bc 22065 genes is too much to plot for the heatmap )
  ## need to double check w/ Andressa what the cutoff should be 
  #DESEQdata.filtered <- DESEQdata[DESEQdata$padj <= 0.01, ]
  
  #Even after trying to filter by 0.01 padj, there was too much data, so I am just going to take the TOP 20 rows for now (20 lowest padj values)
  DESEQdata.filtered <- head(DESEQdata, 20)
  
#SUBSETTING DATA OF INTEREST 
  #1. Keeping only columns of interest (getting rid of statistical columns )
  heatmapdata <- DESEQdata.filtered[, c( 2, 9, 10, 11, 12, 13, 14, 15, 16)]
  #2. Using tidyverse package to reshape df --> long format (compatible w/ heatmap fxn)
  heatmapdata2 <- heatmapdata%>%
    pivot_longer(cols = starts_with("Sample"), 
                 names_to = "Sample", 
                 values_to = "Expression")
  
  ##Set expression values as numeric 
  #heatmapdata2$Expression <- as.numeric(heatmapdata2$Expression)
 
  
  
   
#HEATMAP 
  #creating the plot
  DEheatmap=ggplot(heatmapdata2)+ geom_tile(aes(x=Sample,y=Gene,fill=Expression))+  scale_fill_gradient2(
    #low ="lightblue",
    mid = "lightblue",
    high = "red",
    midpoint = 1943.356,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill")+
    labs(x="Tissue Replicate", y="Gene ID",fill="Normalized Expression (DESeq)") + 
    scale_x_discrete(labels=c("Sample_od1_rep4"="od1",
                              "Sample_od1_rep5"="od1",
                              "Sample_od2_rep2"="od2",
                              "Sample_od2_rep6"="od2",
                              "Sample_od3_rep1"="od3",
                              "Sample_od3_rep5"="od3",
                              "Sample_ovary_rep5"="ovary",
                              "Sample_ovary_rep7"="ovary"))+
    theme(axis.text.x=element_text(angle=20,size=9, vjust=0.88,hjust=0.9),
          axis.text.y=element_text(size=7),
          axis.title.x=element_text(size=11),
          axis.title.y = element_text(size=11),
          plot.margin=margin(0.8,0.8,0.8,0.8,"cm"),legend.position="right")+ ggtitle("Differential Expression Analysis - Reproductive Tissues")
  
  
  
  
  
  
  #printing the plot
  pdf("DEheatmap.pdf", width = 10, height = 10)
  
  print(DEheatmap,
        labels = "Differential Expression Analysis - Reproductive Tissues",
        ncol = 1, nrow = 1, legend = NULL,
        common.legend = TRUE,
        legend.grob = NULL)
  
  
  dev.off()
  

  
  