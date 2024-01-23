##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/pairwise")

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
#install.packages("gridExtra")
library(gridExtra)


#IMPORT DATA 
cntdata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/all_cnt_0625.txt", header=T, na.strings="NA")
phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2.txt", header=T)

#SUBSET phenodata --> 6 pairwise comparisons 
  #SUBSET ALL REPROD. TISSUE SAMPLES (just for easy viewing and checking)
  phenodatarep= subset(phenodata, System == "Reproductive")
  #Subset DFs for each pairwise comparison 
  phenodataod12 <- subset(phenodata, Tissue %in% c("Magnum(od1)", "Isthmus(od2)"))
  phenodataod13 <- subset(phenodata, Tissue %in% c("Magnum(od1)", "Shell_Gland(od3)"))
  phenodataod1o <- subset(phenodata, Tissue %in% c("Magnum(od1)", "Ovary"))
  phenodataod23 <- subset(phenodata, Tissue %in% c("Isthmus(od2)", "Shell_Gland(od3)"))
  phenodataod2o <- subset(phenodata, Tissue %in% c("Isthmus(od2)", "Ovary"))
  phenodataod3o <- subset(phenodata, Tissue %in% c("Shell_Gland(od3)", "Ovary"))

#SUBSET cntdata (gene matrix) --> only samples that matched those from the phenodata of each 6 pairwise comparisons 
  #Isolating only the samples in cntdata [RNAreads] that have matches to the samples listed in each phenodata, isolating only the expressed RNA reads/ tissue pair 
  samplematchesod12=(colnames(cntdata) %in% phenodataod12$Sample_name)
  samplematchesod13=(colnames(cntdata) %in% phenodataod13$Sample_name)
  samplematchesod1o=(colnames(cntdata) %in% phenodataod1o$Sample_name)
  samplematchesod23=(colnames(cntdata) %in% phenodataod23$Sample_name)
  samplematchesod2o=(colnames(cntdata) %in% phenodataod2o$Sample_name)
  samplematchesod3o=(colnames(cntdata) %in% phenodataod3o$Sample_name)
  #create subset gene matrix (cntdata) dfs 
  cntdataod12=cntdata[, samplematchesod12]
  cntdataod13=cntdata[, samplematchesod13]
  cntdataod1o=cntdata[, samplematchesod1o]
  cntdataod23=cntdata[, samplematchesod23]
  cntdataod2o=cntdata[, samplematchesod2o]
  cntdataod3o=cntdata[, samplematchesod3o]
  
  #the isolation will unintentionally get rid of the gene column (col1), so we need to add it back on, but as rownames instead this time . 
  #adding back on the gene column 
    #getting gene names 
      gene = as.data.frame(cntdata$gene)
      colnames(gene) = "gene_id"
    #adding back on the gene names column 
      mRNAdataod12= cbind(gene, cntdataod12)
      mRNAdataod13= cbind(gene, cntdataod13)
      mRNAdataod1o= cbind(gene, cntdataod1o)
      mRNAdataod23= cbind(gene, cntdataod23)
      mRNAdataod2o= cbind(gene, cntdataod2o)
      mRNAdataod3o= cbind(gene, cntdataod3o)
    #converting rownames --> gene_ID values, as opposed to 1,2,3,4,5,etc.. (transform rownames)
    #first getting rid of "gene-" in front of each gene name 
      mRNAdataod12$gene_id <- sub("^gene-", "", mRNAdataod12$gene_id)
      mRNAdataod13$gene_id <- sub("^gene-", "", mRNAdataod13$gene_id)
      mRNAdataod1o$gene_id <- sub("^gene-", "", mRNAdataod1o$gene_id)
      mRNAdataod23$gene_id <- sub("^gene-", "", mRNAdataod23$gene_id)
      mRNAdataod2o$gene_id <- sub("^gene-", "", mRNAdataod2o$gene_id)
      mRNAdataod3o$gene_id <- sub("^gene-", "", mRNAdataod3o$gene_id)
    #transform rownames
      row.names(mRNAdataod12) = mRNAdataod12[,1]
      row.names(mRNAdataod13) = mRNAdataod13[,1]
      row.names(mRNAdataod1o) = mRNAdataod1o[,1]
      row.names(mRNAdataod23) = mRNAdataod23[,1]
      row.names(mRNAdataod2o) = mRNAdataod2o[,1]
      row.names(mRNAdataod3o) = mRNAdataod3o[,1]
    #deleting the extra column ( don't need it anymore since it is the rowname now)
      mRNAdataod12 = mRNAdataod12[,-1]       
      mRNAdataod13 = mRNAdataod13[,-1]  
      mRNAdataod1o = mRNAdataod1o[,-1]  
      mRNAdataod23 = mRNAdataod23[,-1]  
      mRNAdataod2o = mRNAdataod2o[,-1]  
      mRNAdataod3o = mRNAdataod3o[,-1]  
      
#At this point we've created phenodata + mRNAdata for each pairwise comparison 
    
#QUALITY CTRL 
  #FILTER OUT GENES w/ 0 READS 
    #1. calc. rowsums (total # of reads / gene) for mRNAdata 
    od12totalcount=rowSums(mRNAdataod12)
    table(od12totalcount==0)
    
    od13totalcount=rowSums(mRNAdataod13)
    table(od13totalcount==0)
    
    od1ototalcount=rowSums(mRNAdataod1o)
    table(od1ototalcount==0)
    
    od23totalcount=rowSums(mRNAdataod23)
    table(od23totalcount==0)
    
    od2ototalcount=rowSums(mRNAdataod2o)
    table(od2ototalcount==0)
    
    od3ototalcount=rowSums(mRNAdataod3o)
    table(od3ototalcount==0)
  
    #2. create vector of genes that DO NOT have a total rowsume = 0 
    non_zero_genes_od12 <- od12totalcount != 0
    
    non_zero_genes_od13 <- od13totalcount != 0
    
    non_zero_genes_od1o <- od1ototalcount != 0
    
    non_zero_genes_od23 <- od23totalcount != 0
    
    non_zero_genes_od2o <- od2ototalcount != 0
    
    non_zero_genes_od3o <- od3ototalcount != 0
    
    #3. Subset df to filter out 0-read-genes 
    mRNAod12.filtered = mRNAdataod12[non_zero_genes_od12, ]
    
    mRNAod13.filtered = mRNAdataod13[non_zero_genes_od13, ]
    
    mRNAod1o.filtered = mRNAdataod1o[non_zero_genes_od1o, ]
    
    mRNAod23.filtered = mRNAdataod23[non_zero_genes_od23, ]
    
    mRNAod2o.filtered = mRNAdataod2o[non_zero_genes_od2o, ]
    
    mRNAod3o.filtered = mRNAdataod3o[non_zero_genes_od3o, ]

  
  #NORMALIZATION USING TMM [Trimmed-Mean of M-values]
    #1. mRNA.filtered -- [ DGEList ] --> formats as an obj.that is compatible w/ edgeR
    DGEod12 = DGEList(counts = mRNAod12.filtered)
    DGEod13 = DGEList(counts = mRNAod13.filtered)
    DGEod1o = DGEList(counts = mRNAod1o.filtered)
    DGEod23 = DGEList(counts = mRNAod23.filtered)
    DGEod2o = DGEList(counts = mRNAod2o.filtered)
    DGEod3o = DGEList(counts = mRNAod3o.filtered)
    #2. DGE mRNA obj. -- [calcNormFactors...method = TMM]--> calc.s norm. factors via TMM method --> think of it as calc.-ing the # to multiply each library size to bring them to equal proportions [accounting for extreme variation in library size] 
    TMMod12 = calcNormFactors(DGEod12, method = "TMM")
    TMMod13 = calcNormFactors(DGEod13, method = "TMM")
    TMMod1o = calcNormFactors(DGEod1o, method = "TMM")
    TMMod23 = calcNormFactors(DGEod23, method = "TMM")
    TMMod2o = calcNormFactors(DGEod2o, method = "TMM")
    TMMod3o = calcNormFactors(DGEod3o, method = "TMM")
    #3. TMM obj. --[cpm...log=FALSE]--> calc.s Counters Per Million / gene in mRNAdataset (TMM)
      ##the read #s are being converted to CPM units ( standard for DESEQ2)
    cpmod12 = cpm(TMMod12, log = FALSE)
    cpmod13 = cpm(TMMod13, log = FALSE)
    cpmod1o = cpm(TMMod1o, log = FALSE)
    cpmod23 = cpm(TMMod23, log = FALSE)
    cpmod2o = cpm(TMMod2o, log = FALSE)
    cpmod3o = cpm(TMMod3o, log = FALSE)
    #4. cpm obj. --[log(cpm+1, base=2)]--> +1 to all CPMs (stabilizes data) + calc.s base-2 log of each result --> obj. w/ log-transformed (TMM-normalized) CPM values used for DE analysis
      #logging the values will fix the distribution of the data, so that it is easier for softwares to work with. 
        #note: the best method to use for normalization is always based off the distribution type of your data, it is now a method that you pick just bc you like it.
          #+1 gets rid of (-) #s, we usually stick to a base = 2 OR 10
    log.cpmod12 = log(cpmod12 + 1, base = 2)
    log.cpmod13 = log(cpmod13 + 1, base = 2)
    log.cpmod1o = log(cpmod1o + 1, base = 2)
    log.cpmod23 = log(cpmod23 + 1, base = 2)
    log.cpmod2o = log(cpmod2o + 1, base = 2)
    log.cpmod3o = log(cpmod3o + 1, base = 2)

#PCA [CHOSE SVD PCA - want to look at important patterns in data, instead of choosing spectral decomp. (looking @ relationship b/w genes themselves)]
##"PCA also tells you which variables (genes, in this case) contribute most to the variability in your data. This can help identify genes that are most relevant to the differences between tissues."
    
    pdf("PCA.pdf", width=10,height=7)
    
    par(mfrow = c(3, 2))
    
    nb.cols <- 20
    mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
    
    pca.countsod12 = prcomp(t(log.cpmod12), scale=F)
    fviz_pca_ind(pca.countsod12,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod12$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    pca.countsod13 = prcomp(t(log.cpmod13), scale=F)
    fviz_pca_ind(pca.countsod13,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod13$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    pca.countsod1o = prcomp(t(log.cpmod1o), scale=F)
    fviz_pca_ind(pca.countsod1o,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod1o$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    pca.countsod23 = prcomp(t(log.cpmod23), scale=F)
    fviz_pca_ind(pca.countsod23,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod23$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    pca.countsod2o = prcomp(t(log.cpmod2o), scale=F)
    fviz_pca_ind(pca.countsod2o,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod2o$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    pca.countsod3o = prcomp(t(log.cpmod3o), scale=F)
    fviz_pca_ind(pca.countsod3o,
                 pointsize=2,
                 geom.ind = ("point"), # show points only (nbut not "text")
                 habillage = phenodataod3o$Tissue, # color by groups
                 palette = mycolors,
                 legend.title = "Tissue",
                 mean.point = FALSE
    )
    
    dev.off()
    
    
#BOXPLOT (showing expression distribution post-TMM Normalization)
#Code Template: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/

    #1 Setting the new Color Palette for the boxplots 
    nsamplesod12 <- ncol(log.cpmod12)
    colod12 <- brewer.pal(nsamplesod12, "Paired")
    nsamplesod13 <- ncol(log.cpmod13)
    colod13 <- brewer.pal(nsamplesod13, "Paired")
    nsamplesod1o <- ncol(log.cpmod1o)
    colod1o <- brewer.pal(nsamplesod1o, "Paired")
    nsamplesod23 <- ncol(log.cpmod23)
    colod23 <- brewer.pal(nsamplesod23, "Paired")
    nsamplesod2o <- ncol(log.cpmod2o)
    colod2o <- brewer.pal(nsamplesod2o, "Paired")
    nsamplesod3o <- ncol(log.cpmod3o)
    colod3o <- brewer.pal(nsamplesod3o, "Paired")
    
    #2 making the boxplots 
    
    pdf("ExpressionDistr.pdf", width=17,height=10)
    
    par(mfrow = c(2, 3))
    
    distrod12=boxplot(log.cpmod12, las=0, col=colod12, main="OD12 Expression Distribution")
    title(main="", ylab="Log-cpm", xlab = "Samples")
    
    distrod13=boxplot(log.cpmod13, las=0, col=colod13, main="OD13 Expression Distribution")
    title(main="", ylab="Log-cpm")
    
    distrod1o=boxplot(log.cpmod1o, las=0, col=colod1o, main="OD1O Expression Distribution")
    title(main="", ylab="Log-cpm")
    
    distrod23=boxplot(log.cpmod23, las=0, col=colod23, main="OD23 Expression Distribution")
    title(main="", ylab="Log-cpm")
    
    distrod2o=boxplot(log.cpmod2o, las=0, col=colod2o, main="OD2O Expression Distribution")
    title(main="", ylab="Log-cpm")
    
    distrod3o=boxplot(log.cpmod3o, las=0, col=colod3o, main="OD3O Expression Distribution")
    title(main="", ylab="Log-cpm")
    
    dev.off()
    
  
   
   