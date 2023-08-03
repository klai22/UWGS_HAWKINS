#setwd/load & install newlibraries 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/GlobalDMR")
#install.packages("pheatmap")
#install.packages("svglite")
#install.packages("ggplotify")
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(svglite)
library(viridis)
library(corrplot)
library(gridExtra)
library(ggplotify)
library(ggpubr)

###HEATMAPs
##od12
#Import / Edit Data 
od12 = read.csv("od1-od2_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od12sig=head(od12,n=1000)
head(od12sig)
#creating new df isolating only columns of interest (methylation lvls)
od12met= od12sig[,c("mC_A","mC_B")]
colnames(od12met)= c("magnum", "isthmus")
head(od12met)

pheatmap(od12met)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od12global=pheatmap(od12met, heat.colors, cluster_rows = T, show_rownames=F,
         border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("Low", "Medium", "High"),
         main = "Magnum vs. Isthmus DMRs",
         cellwidth = 100
)

##od13
#Import / Edit Data 
od13 = read.csv("od1-od3_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od13sig=head(od13,n=1000)
head(od13sig)
#creating new df isolating only columns of interest (methylation lvls)
od13met= od13sig[,c("mC_A","mC_B")]
colnames(od13met)= c("magnum", "shell gland")
head(od13met)

pheatmap(od13met)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od13global=pheatmap(od13met, heat.colors, cluster_rows = T, show_rownames=F,
                    border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
                    legend_breaks = c(-2, 0, 2),
                    legend_labels = c("Low", "Medium", "High"),
                    main = "Magnum vs. Shell Gland DMRs",
                    cellwidth = 100
)

##od1o
#Import / Edit Data 
od1o = read.csv("od1-ovary_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od1osig=head(od1o,n=1000)
head(od1osig)
#creating new df isolating only columns of interest (methylation lvls)
od1omet= od1osig[,c("mC_A","mC_B")]
colnames(od1omet)= c("magnum", "ovary")
head(od1omet)

pheatmap(od1omet)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od1oglobal=pheatmap(od1omet, heat.colors, cluster_rows = T, show_rownames=F,
                    border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
                    legend_breaks = c(-2, 0, 2),
                    legend_labels = c("Low", "Medium", "High"),
                    main = "Magnum vs. Ovary DMRs",
                    cellwidth = 100
)

##od23
#Import / Edit Data 
od23 = read.csv("od2-od3_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od23sig=head(od23,n=1000)
head(od23sig)
#creating new df isolating only columns of interest (methylation lvls)
od23met= od23sig[,c("mC_A","mC_B")]
colnames(od23met)= c("isthmus", "shell gland")
head(od23met)

pheatmap(od23met)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od23global=pheatmap(od23met, heat.colors, cluster_rows = T, show_rownames=F,
                    border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
                    legend_breaks = c(-2, 0, 2),
                    legend_labels = c("Low", "Medium", "High"),
                    main = "Isthmus vs. Shell gland DMRs",
                    cellwidth = 100
)

##od2o
#Import / Edit Data 
od2o = read.csv("od2-ovary_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od2osig=head(od2o,n=1000)
head(od2osig)
#creating new df isolating only columns of interest (methylation lvls)
od2omet= od2osig[,c("mC_A","mC_B")]
colnames(od2omet)= c("isthmus", "ovary")
head(od2omet)

pheatmap(od2omet)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od2oglobal=pheatmap(od2omet, heat.colors, cluster_rows = T, show_rownames=F,
                    border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
                    legend_breaks = c(-2, 0, 2),
                    legend_labels = c("Low", "Medium", "High"),
                    main = "Isthmus vs. Ovary DMRs",
                    cellwidth = 100
)

##od3o
#Import / Edit Data 
od3o = read.csv("od3-ovary_filter_DMR_mod.csv", header=T, row.names = 1)
#isolating top 1000 DMRs (threshold)
od3osig=head(od3o,n=1000)
head(od3osig)
#creating new df isolating only columns of interest (methylation lvls)
od3omet= od3osig[,c("mC_A","mC_B")]
colnames(od3omet)= c("shell gland", "ovary")
head(od3omet)

pheatmap(od3omet)

#Plot prep & creation 
#sample_col= data.frame(group=data[,'group'])
heat.colors <- brewer.pal(8, "Reds")
od3oglobal=pheatmap(od3omet, heat.colors, cluster_rows = T, show_rownames=F,
                    border_color=NA, fontsize = 8, legend = TRUE,scale = "column",
                    legend_breaks = c(-2, 0, 2),
                    legend_labels = c("Low", "Medium", "High"),
                    main = "Shell Gland vs. Ovary DMRs",
                    cellwidth = 100
)

# convert pheatmap objects to ggplot objects
od12globalgrob <- as.ggplot(od12global)
od13globalgrob <- as.ggplot(od13global)
od1oglobalgrob <- as.ggplot(od1oglobal)
od23globalgrob <- as.ggplot(od23global)
od2oglobalgrob <- as.ggplot(od2oglobal)
od3oglobalgrob <- as.ggplot(od3oglobal)

# save to pdf
pdf("globalDMRs.pdf", width=15,height=10)
print(ggarrange(od12globalgrob,od13globalgrob,od1oglobalgrob,od23globalgrob,od2oglobalgrob,od3oglobalgrob,
             ncol=3,nrow=2))
dev.off()


###CORRELATION PLOT 
data1=read.table("data.c_read_Count-filter2.txt", header=T)
data1[1:4,1:4]
colnames(data1)=c("chr.base","od1","od2","od3"," Ovary")
data1[1:4,1:4]
## row.names
row.names(data1) = data1[,1]
data1 = data1[,-1]
data1[1:4,1:4]
##filter zero
#remove genes non-expressed and low expressed
## total counts per gene
Totalcounts = rowSums(data1)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(data1)==0
data2 = data1[!rm,]
dim(data2)#
##log.fpkm
log.data=log(data2 +1 , base=2)
dim(log.data)
heat.colors <- brewer.pal(11, "RdBu")
#heat.colors = viridis(250)
pheatmap(log.data, heat.colors, cluster_rows = T, cluster_cols = T,show_rownames=F,
         border_color=NA, fontsize = 8, legend = TRUE,scale = "column", method = "complete", kmeans_k = NA,
         legend_breaks = c(-2, 0, 2), clustering_distance_rows = "euclidean",  clustering_distance_cols = "euclidean",
         legend_labels = c("Low", "Medium", "High"))
corr_mat=cor(data2, method = "s")
correlationplot=corrplot(corr_mat, method = 'number', type="lower") # colorful number

# save to pdf
pdf("DMRCorrelationPlot.pdf", width=15,height=10)
print(corrplot(corr_mat, method = 'number', type="lower"))
dev.off()



