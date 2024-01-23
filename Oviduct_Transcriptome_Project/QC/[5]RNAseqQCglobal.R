#Kenneth Lai 
##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/global")

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
#install.packages("rsvg")
library(rsvg)
library(grDevices)
library(MASS)



#IMPORT DATA 
cntdata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/all_cnt_0625.txt", header=T, na.strings="NA")
phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2.txt", header=T)

#ISOLATING EXPRESSED RNA 
##Isolating only the samples in cntdata [RNAreads] that have matches to the samples listed in phenodata, isolating only the expressed RNA reads
samplematches=(colnames(cntdata) %in% phenodata$Sample_name)
#table(samplematches)
mRNAdata= cntdata[, samplematches]
#the isolation will unintentionally get rid of the gene column (col1), so we need to add it back on, but as rownames instead this time . 
#adding back on the gene column 
gene = as.data.frame(cntdata$gene)
colnames(gene) = "gene_id"
mRNAdata= cbind(gene, mRNAdata)
#converting rownames --> gene_ID values, as opposed to 1,2,3,4,5,etc.. (transform rownames)
#first getting rid of "gene-" in front of each gene name 
mRNAdata$gene_id <- sub("^gene-", "", mRNAdata$gene_id)
#transform rownames
row.names(mRNAdata) = mRNAdata[,1]
#deleting the extra column ( don't need it anymore since it is the rowname now)
mRNAdata = mRNAdata[,-1]          

#CREATING SUBDFS for our tissues of interest 
#subsetting phenodata/tissue 
phenodatareproductive= subset(phenodata, System == "Reproductive")

#subsetting mRNAdata/tissue 
# Get the column names (sampleID) / tissue type 
reproductivesamples <- phenodatareproductive$Sample_name

# Subset DF2 based on the column names
mRNAreproductive <- mRNAdata[, names(mRNAdata) %in% reproductivesamples]
  

#At this point we've created phenodata + mRNAdata for reproductive tissues only  
    
#QUALITY CTRL 
  #FILTER OUT GENES w/ 0 READS 
    #1. calc. rowsums (total # of reads / gene) for mRNAdata 
    totalcounts=rowSums(mRNAreproductive)
    table(totalcounts==0)
  
    #2. create vector of genes that DO NOT have a total rowsum = 0 
    non_zero_genes <- totalcounts != 0
    
    #3. Subset df to filter out 0-read-genes 
    mRNAreproductive.filtered = mRNAreproductive[non_zero_genes, ]
    
  
  
  #NORMALIZATION USING TMM [Trimmed-Mean of M-values]
    #1. mRNA.filtered -- [ DGEList ] --> formats as an obj.that is compatible w/ edgeR
    DGEreproductive = DGEList(counts = mRNAreproductive.filtered)
    #2. DGE mRNA obj. -- [calcNormFactors...method = TMM]--> calc.s norm. factors via TMM method --> think of it as calc.-ing the # to multiply each library size to bring them to equal proportions [accounting for extreme variation in library size] 
    TMMreproductive = calcNormFactors(DGEreproductive, method = "TMM")
    #3. TMM obj. --[cpm...log=FALSE]--> calc.s Counters Per Million / gene in mRNAdataset (TMM)
      ##the read #s are being converted to CPM units ( standard for DESEQ2)
    cpmreproductive = cpm(TMMreproductive, log = FALSE)
    #4. cpm obj. --[log(cpm+1, base=2)]--> +1 to all CPMs (stabilizes data) + calc.s base-2 log of each result --> obj. w/ log-transformed (TMM-normalized) CPM values used for DE analysis
      #logging the values will fix the distribution of the data, so that it is easier for softwares to work with. 
        #note: the best method to use for normalization is always based off the distribution type of your data, it is now a method that you pick just bc you like it.
          #+1 gets rid of (-) #s, we usually stick to a base = 2 OR 10
    log.cpmreproductive = log(cpmreproductive + 1, base = 2)

#PCA [CHOSE SVD PCA - want to look at important patterns in data, instead of choosing spectral decomp. (looking @ relationship b/w genes themselves)]
##"PCA also tells you which variables (genes, in this case) contribute most to the variability in your data. This can help identify genes that are most relevant to the differences between tissues."
  #PCA-alltissues 
    pdf("PCA_reproductive_tissues.pdf", width = 8, height = 6)
    par(mar = c(8, 4, 4, 4))
    par(cex = 0.8)
    
    # Define the number of colors you want
    nb.cols <- 20
    mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
    pca.counts2 = prcomp(t(log.cpmreproductive), scale=F)
    fviz_pca_ind(pca.counts2,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodatareproductive$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    dev.off() 
    
#MDS[Multi-Dimensional Scaling] (we didn't like how PC1 clustered the tissues together as much, so we will try using MDS instead of PCA to see if results align better w/ dendogram results)

    #1. Apply limma-plotMDS to do actual MDS analysis 
    mds <- data.frame(plotMDS(cpmreproductive)[c("x", "y")])
    #2. Merge (cbind) mds + pheno data 
    mds <- cbind(mds, phenodatareproductive)
    #3. Creating MDS as a plot 
    mdsplot <- ggplot(mds) +
      aes(x=x, y=y, color=Tissue, shape=Tissue) +
      xlab("PC1") + ylab("PC2") +
      geom_point(size=3) +
      coord_fixed(ratio=1) +
      ggtitle("MDS Plot of Gene expression by Tissue") +
      scale_color_manual(values=c("blueviolet", "yellow",  "green","blue")) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(fill= "transparent"))
    #4. Save plot as pdf 
    pdf("MDS_reproductive_tissues.pdf", width = 8, height = 6)
    par(mar = c(8, 4, 4, 4))
    par(cex = 0.8)
    mdsplot
    dev.off()
    
#BOXPLOT (showing expression distribution post-TMM Normalization)
#Code Template: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/

  #A. POST NORMALIZATION PLOT 
    
    #1 Setting the new Color Palette for the boxplots 
    nsamples <- ncol(log.cpmreproductive)
    col <- brewer.pal(nsamples, "Paired")

    #2 making the boxplots 
    
    pdf("ExpressionDistr_reproductive_tissues.pdf", width=17,height=10)
    
    distrreproductive=boxplot(log.cpmreproductive, las=0, col=col, main="Expression Distribution")
    title(main="", ylab="Log-cpm", xlab = "Samples")
    
    dev.off()
    
  #B. PRE NORMALIZATION PLOT 
    nsamples2 <- ncol(mRNAreproductive.filtered)
    col2 <- brewer.pal(nsamples2, "Paired")
    
    #2 making the boxplots 
    
    pdf("ExpressionDistr_reproductive_tissues[pre].pdf", width=17,height=10)
    
    distrreproductive2=boxplot(mRNAreproductive.filtered, las=0, col2=col2, main="Expression Distribution")
    title(main="", ylab="reads", xlab = "Samples")
    
    dev.off()
 
#DENDROGRAM (clustering samples by similarity of gene expression)
##taken from Andressa's template code 
  #1 converting to df + transposing (switch rows/columns place) to make compatible w/ hclust software
    datExpr = t(as.data.frame(log.cpmreproductive))
    
  #2 creating pdf 
    svglite("dendogram_reproducitve_tissues.svg", width=10, height=8, pointsize = 8)
    par(mar = c(8,4,4,4))
    par (cex=1.0)
  #3 hclust (hierarchical clustering) datExpr using the "average" linkage method 
    sampleTree = hclust(dist(datExpr), method = "average")
  #4 plotting result 
    plot(sampleTree, main = "Cluster Dendrogram", sub="Reproductive Tissues", 
         xlab="Replicate Samples", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2.5)
    dev.off()
  #5 converting svg file --[rsvg R package]--> pdf file
    rsvg::rsvg_png("dendogram_reproducitve_tissues.svg", file = "dendogram_reproducitve_tissues.png")
    
    png("dendrogram_reproductive_tissues.pdf", width = 10, height = 8, units = "in", res = 300)
    rasterImage(readPNG("dendogram_reproducitve_tissues.png"), 0, 0, 10, 8)
    dev.off()
    
    
    