#Kenneth Lai 

#setwd
setwd("/Users/kennethlai/desktop/HAWKINS/FAANG/RNAseq/Volcano")
#libraries 
library(ggplot2)
library("ggpubr")
library(RColorBrewer)
library(pheatmap)
library(svglite)
library(viridis)
library(corrplot)
library(gridExtra)
library(ggplotify)
library(ggpubr)
library("wesanderson")
library("ComplexHeatmap")
library("dplyr")
library("tibble")
library(volcanoPlot)

#Import Data (DESEQ RESULTS)
od12 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od12]DESEQ-results7.csv")
od13 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od13]DESEQ-results7.csv")
od1o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od1o]DESEQ-results7.csv")
od23 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od23]DESEQ-results7.csv")
od2o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od2o]DESEQ-results7.csv")
od3o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od3o]DESEQ-results7.csv")
global = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[global]DESEQ-results7.csv")

#Adding a -log10(p-value) column for each df 
od12$neglog10pvalues <- -log10(od12$pvalue)
od13$neglog10pvalues <- -log10(od13$pvalue)
od1o$neglog10pvalues <- -log10(od1o$pvalue)
od23$neglog10pvalues <- -log10(od23$pvalue)
od2o$neglog10pvalues <- -log10(od2o$pvalue)
od3o$neglog10pvalues <- -log10(od3o$pvalue)



#PLOTTING VOLCANO PLOTS 
  #identifying top -log10(p-values)
    top_positiveod12 <- od12[od12$log2FoldChange > 0, ][order(-od12$neglog10pvalues), ][1:10, ]
    top_negativeod12 <- od12[od12$log2FoldChange < 0, ][order(-od12$neglog10pvalues), ][1:10, ]
    top_positiveod13 <- od13[od13$log2FoldChange > 0, ][order(-od13$neglog10pvalues), ][1:10, ]
    top_negativeod13 <- od13[od13$log2FoldChange < 0, ][order(-od13$neglog10pvalues), ][1:10, ]
    top_positiveod1o <- od1o[od1o$log2FoldChange > 0, ][order(-od1o$neglog10pvalues), ][1:10, ]
    top_negativeod1o <- od1o[od1o$log2FoldChange < 0, ][order(-od1o$neglog10pvalues), ][1:10, ]
    top_positiveod23 <- od23[od23$log2FoldChange > 0, ][order(-od23$neglog10pvalues), ][1:10, ]
    top_negativeod23 <- od23[od23$log2FoldChange < 0, ][order(-od23$neglog10pvalues), ][1:10, ]
    top_positiveod2o <- od2o[od2o$log2FoldChange > 0, ][order(-od2o$neglog10pvalues), ][1:10, ]
    top_negativeod2o <- od2o[od2o$log2FoldChange < 0, ][order(-od2o$neglog10pvalues), ][1:10, ]
    top_positiveod3o <- od3o[od3o$log2FoldChange > 0, ][order(-od3o$neglog10pvalues), ][1:10, ]
    top_negativeod3o <- od3o[od3o$log2FoldChange < 0, ][order(-od3o$neglog10pvalues), ][1:10, ]


  #plotting 
    #RE-COLORED, I think ordering got mixed up (couldn't tell which tissue was positive or negative per pair, I reassessed to match the DE heatmap)
    #Color Code: ("Magnum" = "blue", "Isthmus" = "yellow","Shell_Gland" = "aquamarine2","Ovary" = "palevioletred2")
    od12volcano=ggplot(od12, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#F87008", "#27DAD8"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Magnum vs. Isthmus Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod12, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod12, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)
    
    od13volcano=ggplot(od13, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#F87008", "#22E11D"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Magnum vs. Shell Gland Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod13, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod13, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)
    
    od1ovolcano=ggplot(od1o, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#F87008", "#F215EE"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Magnum vs. Ovary Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod1o, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod1o, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)
    
    od23volcano=ggplot(od23, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#22E11D", "#27DAD8"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Shell Gland vs. Isthmus Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod23, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod23, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)
    
    od2ovolcano=ggplot(od2o, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#F215EE", "#27DAD8"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Ovary vs. Ithmus Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod2o, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod2o, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)
    
    od3ovolcano=ggplot(od3o, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = ifelse(pvalue > 0.05, "grey",
                                    ifelse(log2FoldChange > 0, "#22E11D", "#F215EE"))), size = 2) +
      scale_color_identity() +
      labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Shell Gland vs. Ovary Differential Expression") +
      theme_minimal() +
      # Add labels for the top 10 positive log2FC values
      geom_text(data = top_positiveod3o, aes(label = Gene), vjust = -0.5, hjust = 1, color = "black", size = 3) +
      # Add labels for the top 10 negative log2FC values
      geom_text(data = top_negativeod3o, aes(label = Gene), vjust = -0.5, hjust = -0.5, color = "black", size = 3)


# save to pdf
    #making plot list 
    volcanoplots=list(od12volcano,od13volcano,od1ovolcano,od23volcano,od2ovolcano,od3ovolcano)

    #making pdf 
    pdf("[5]Volcano_Plots.pdf", width = 25, height = 13)
    ggarrange(plotlist = volcanoplots, ncol = 3, nrow = 2)
    dev.off()
    