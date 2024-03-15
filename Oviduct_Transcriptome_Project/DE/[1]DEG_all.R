#Kenneth Lai 

#setwd
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/DEFigures2")
#libraries 
library(ggplot2)
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
library(devtools)
library(scales)
source("/Users/kennethlai/desktop/HAWKINS/FAANG/RNAseq/colorpalette/RHodor_colors/R/RHodor_palettes.R")

#Import Data (DESEQ RESULTS)
od12 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od12]DESEQ-results7.csv")
od13 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od13]DESEQ-results7.csv")
od1o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od1o]DESEQ-results7.csv")
od23 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od23]DESEQ-results7.csv")
od2o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od2o]DESEQ-results7.csv")
od3o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[od3o]DESEQ-results7.csv")
global = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE/[global]DESEQ-results7.csv")

#CLEANING DATA (taking only most signifigant genes)
#Including all data this time
  #order rows (greatest-->smallest) based on |logFC|
od12 <- od12 %>% arrange(desc(abs(log2FoldChange)))
od13 <- od13 %>% arrange(desc(abs(log2FoldChange)))
od1o <- od1o %>% arrange(desc(abs(log2FoldChange)))
od23 <- od23 %>% arrange(desc(abs(log2FoldChange)))
od2o <- od2o %>% arrange(desc(abs(log2FoldChange)))
od3o <- od3o %>% arrange(desc(abs(log2FoldChange)))
  #Select the top 20 rows (20 highest |logFC|)
#od12 <- head(od12, 20)
#od13 <- head(od13, 20)
#od1o <- head(od1o, 20)
#od23 <- head(od23, 20)
#od2o <- head(od2o, 20)
#od3o <- head(od3o, 20)

#Creating new dfs w/ only columns of interest (expression data + genes ) p.1
od12_filtered= od12[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]
od13_filtered= od13[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]
od1o_filtered= od1o[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]
od23_filtered= od23[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]
od2o_filtered= od2o[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]
od3o_filtered= od3o[c("Gene","Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2","Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5","Sample_ovary_rep5","Sample_ovary_rep7")]

#Editing dfs w/ new columns that AVERAGE each of the replicate pairs' expression values 
column_pair1 <- c("Sample_od1_rep4","Sample_od1_rep5")
column_pair2 <- c("Sample_od2_rep2","Sample_od2_rep6")
column_pair3 <- c("Sample_od3_rep1","Sample_od3_rep5")
column_pair4 <- c("Sample_ovary_rep5","Sample_ovary_rep7")
              
od12_filtered$od1 <- rowMeans(od12_filtered[column_pair1])
od12_filtered$od2 <- rowMeans(od12_filtered[column_pair2])
od12_filtered$od3 <- rowMeans(od12_filtered[column_pair3])
od12_filtered$ovary <- rowMeans(od12_filtered[column_pair4])
od13_filtered$od1 <- rowMeans(od13_filtered[column_pair1])
od13_filtered$od2 <- rowMeans(od13_filtered[column_pair2])
od13_filtered$od3 <- rowMeans(od13_filtered[column_pair3])
od13_filtered$ovary <- rowMeans(od13_filtered[column_pair4])
od1o_filtered$od1 <- rowMeans(od1o_filtered[column_pair1])
od1o_filtered$od2 <- rowMeans(od1o_filtered[column_pair2])
od1o_filtered$od3 <- rowMeans(od1o_filtered[column_pair3])
od1o_filtered$ovary <- rowMeans(od1o_filtered[column_pair4])
od23_filtered$od1 <- rowMeans(od23_filtered[column_pair1])
od23_filtered$od2 <- rowMeans(od23_filtered[column_pair2])
od23_filtered$od3 <- rowMeans(od23_filtered[column_pair3])
od23_filtered$ovary <- rowMeans(od23_filtered[column_pair4])
od2o_filtered$od1 <- rowMeans(od2o_filtered[column_pair1])
od2o_filtered$od2 <- rowMeans(od2o_filtered[column_pair2])
od2o_filtered$od3 <- rowMeans(od2o_filtered[column_pair3])
od2o_filtered$ovary <- rowMeans(od2o_filtered[column_pair4])
od3o_filtered$od1 <- rowMeans(od3o_filtered[column_pair1])
od3o_filtered$od2 <- rowMeans(od3o_filtered[column_pair2])
od3o_filtered$od3 <- rowMeans(od3o_filtered[column_pair3])
od3o_filtered$ovary <- rowMeans(od3o_filtered[column_pair4])

#Creating new dfs w/ only columns of interest (AVEREAGED expression data (b/w replicates) + genes ) p.2
od12_avg= od12_filtered[c("Gene","od1","od2","od3","ovary")]
od13_avg= od13_filtered[c("Gene","od1","od2","od3","ovary")]
od1o_avg= od1o_filtered[c("Gene","od1","od2","od3","ovary")]
od23_avg= od23_filtered[c("Gene","od1","od2","od3","ovary")]
od2o_avg= od2o_filtered[c("Gene","od1","od2","od3","ovary")]
od3o_avg= od3o_filtered[c("Gene","od1","od2","od3","ovary")]

#Merging --> master df + eliminating duplicate genes 
masterDEG <- bind_rows(od12_avg,od13_avg,od1o_avg,od23_avg,od2o_avg,od3o_avg) %>%
  distinct(Gene, .keep_all = TRUE)

#Renaming Sample / Tissue Names 
colnames(masterDEG)[colnames(masterDEG) == "od1"] <- "Magnum(od1)"
colnames(masterDEG)[colnames(masterDEG) == "od2"] <- "Isthmus(od2)"
colnames(masterDEG)[colnames(masterDEG) == "od3"] <- "Shell Gland(od3)"
colnames(masterDEG)[colnames(masterDEG) == "ovary"] <- "Ovary"

#making compatible w/ heat map 
rownames(masterDEG)=masterDEG$Gene
masterDEG=masterDEG[,-1]

#PLOTTING HEATMAPS 
  #(Andressa) requested that row names be omitted (too many)

colors <- hodor_pal("soho") (10)
samplesnames <- c("Magnum(od1)", "Isthmus(od2)", "Shell Gland(od3)", "Ovary")

DEG=Heatmap(
  masterDEG,
  show_column_names = TRUE,
  column_title = "Differential Expression",
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  name = "Expression Level (DESEQ2)",
  row_title = "Signifigant Differentially Expressed Genes",
  row_names_gp = gpar(fontsize = 6, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  col = colors,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  ),
  row_names_side = "left",
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_order = samplesnames)


# save to pdf
pdf("[4]DEG_all.pdf", width=10,height=8)
print(DEG)
dev.off()
