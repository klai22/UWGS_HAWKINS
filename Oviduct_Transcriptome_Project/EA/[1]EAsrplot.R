#Kenneth Lai
# setwd 
setwd("/Users/kennethlai/Desktop/HAWKINS/FAANG/RNAseq/EA/6(SRPlot)")

#load packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("biomaRt")
library("tidyr")


#import data 
od12 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od12]DESEQ-results7.csv")
od13 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od13]DESEQ-results7.csv")
od1o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od1o]DESEQ-results7.csv")
od23 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od23]DESEQ-results7.csv")
od2o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od2o]DESEQ-results7.csv")
od3o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od3o]DESEQ-results7.csv")

#Isolating only columns for gene list + logFC (that's SRPlot's input )
od12=od12[,c("Gene","log2FoldChange")]
od13=od13[,c("Gene","log2FoldChange")]
od1o=od1o[,c("Gene","log2FoldChange")]
od23=od23[,c("Gene","log2FoldChange")]
od2o=od2o[,c("Gene","log2FoldChange")]
od3o=od3o[,c("Gene","log2FoldChange")]
  

#Saving list 
    #write.csv(od12, "od12genelist.csv", row.names = FALSE)
    #write.csv(od13, "od13genelist.csv", row.names = FALSE)
    #write.csv(od1o, "od1ogenelist.csv", row.names = FALSE)
    #write.csv(od23, "od23genelist.csv", row.names = FALSE)
    #write.csv(od2o, "od2ogenelist.csv", row.names = FALSE)
    #write.csv(od3o, "od3ogenelist.csv", row.names = FALSE)
    
#Saving as point-separated values (comaptible w/ SRplot)
    write.table(od12, file = "od12genelist.psv", sep = ".", row.names = FALSE)
    write.table(od13, file = "od13genelist.psv", sep = ".", row.names = FALSE)
    write.table(od1o, file = "od1ogenelist.psv", sep = ".", row.names = FALSE)
    write.table(od23, file = "od23genelist.psv", sep = ".", row.names = FALSE)
    write.table(od2o, file = "od2ogenelist.psv", sep = ".", row.names = FALSE)
    write.table(od3o, file = "od3ogenelist.psv", sep = ".", row.names = FALSE)
    