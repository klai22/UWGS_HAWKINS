#Kenneth Lai
#setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/DE")

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
cntdata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/all_cnt_0625.txt", header=T, na.strings="NA")
phenodata= read.table("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/QC/data/pheno_2(edit).txt", header=T)

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
  #Get the column names (sampleID) / tissue type 
    reproductivesamples <- phenodatareproductive$Sample_name

  #Subset DF2 based on the column names
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

#DE (Differential Expression w/ DESEQ2)
  #Defining the model / selecting inputs dfs 
    #1. set phenodatareproductive$Tissue as a factor (categorial variable) when doing DE
      phenodatareproductive$Tissue = as.factor(phenodatareproductive$Tissue)
    
    #2. Merging mRNAreproductive.filtered + phenodatareproductive dfs --> DESeqDataSet format (1 data structure compatible w/ DESEQ2)
      ##design = ~Tissue specifies that we want to compare DE b/w Tissues (group samples by $Tissue)
      #dds = DeseqDataSet
    ddsreproductive<-DESeqDataSetFromMatrix(countData=mRNAreproductive.filtered, colData= phenodatareproductive, design= ~ Tissue)
    
    #3. Extracting all $Tissue variables in our new dds df (checking)
    colData(ddsreproductive)$Tissue
    
  #DIFFERENTIAL EXPRESSION ANALYSIS (applying DESEQ2 )
    #1. Defining $Tissue (telling DESEQ2) as the categorical-variable as a factor of interest
    design(ddsreproductive) <- formula(~ Tissue)
    #2a. Applying DESEQ2 to do the actual DE, disabling beta prior method
    ddsreproductive <- DESeq(ddsreproductive, betaPrior=FALSE)
    #2b. Checking names of pairwise comparisons done (default: "Variable_Level_vs_Reference_Level,". ref = od2 bc alphabetical order)
    resultsNames(ddsreproductive)
    
    #3. Retrieves results/pairwisecomparison --> apply FDR (false-discovery-rate) method to adjust p-values (already calc.-ed in DESeq step) to help ID true DE genes vs. false (+)s 
    resreproductive<-results(ddsreproductive,pAdjustMethod = "fdr" ) #global 
    resod12<-results(ddsreproductive,name="Tissue_Magnum.od1_vs_Isthmus.od2", pAdjustMethod = "fdr" ) 
    #resod13<-results(ddsreproductive,name="Tissue_Magnum.od1_vs_Shell_Gland.od3", pAdjustMethod = "fdr" ) #MISSING from contrasts [Error: subscript contains invalid names]
    #resod1o<-results(ddsreproductive,name="Tissue_Magnum.od1_vs_Ovary", pAdjustMethod = "fdr" ) #MISSING from contrasts [Error: subscript contains invalid names]
    resod23<-results(ddsreproductive,name="Tissue_Shell_Gland.od3_vs_Isthmus.od2", pAdjustMethod = "fdr" )
    resod2o<-results(ddsreproductive,name="Tissue_Ovary_vs_Isthmus.od2", pAdjustMethod = "fdr" )
    #resod3o<-results(ddsreproductive,name="Tissue_Shell_Gland.od3_vs_Ovary", pAdjustMethod = "fdr" ) #MISSING from contrasts [Error: subscript contains invalid names]
    
    #adding into variable_lvl vs. variable_lvl missing contrasts [od13,od1o,od3o] --> https://support.bioconductor.org/p/111685/
    resod13<-results(ddsreproductive,contrast = c("Tissue","Magnum.od1","Shell_Gland.od3"), pAdjustMethod = "fdr" ) 
    resod1o<-results(ddsreproductive,contrast = c("Tissue","Magnum.od1","Ovary"), pAdjustMethod = "fdr" ) 
    resod3o<-results(ddsreproductive,contrast = c("Tissue","Shell_Gland.od3","Ovary"), pAdjustMethod = "fdr" ) 
  
    resod12test<-results(ddsreproductive,contrast = c("Tissue","Magnum.od1","Isthmus.od2"), pAdjustMethod = "fdr" ) #testing that the code I used for variable lvl comparisons yield same results as refernece lvl comparisons
    #4.Extract metadata associated from dds results, asking to include column names (TRUE)
    #mcols(resod12, use.names=TRUE)
    #5. Reorder results (dds) in order of increasing adjp (adjusted-p-values). Lower padj = more signifigant DE for said gene! 
    resreproductive<-resreproductive[order(resreproductive$padj),]
    resod12<-resod12[order(resod12$padj),]
    resod13<-resod13[order(resod13$padj),]
    resod1o<-resod1o[order(resod1o$padj),]
    resod23<-resod23[order(resod23$padj),]
    resod2o<-resod2o[order(resod2o$padj),]
    resod3o<-resod3o[order(resod3o$padj),]
    
    resod12test<-resod12test[order(resod12test$padj),]
    ##print out first 21 lines of resulting ordered dds 
    #head(resod12, 21)
    #6. Select only rows wheree log fold change >= abs.value of 1--> (logFC >=|1|)
    resreproductive <- subset(resreproductive, abs(log2FoldChange) >=1)
    resod12 <- subset(resod12, abs(log2FoldChange) >=1)
    resod13 <- subset(resod13, abs(log2FoldChange) >=1)
    resod1o <- subset(resod1o, abs(log2FoldChange) >=1)
    resod23 <- subset(resod23, abs(log2FoldChange) >=1)
    resod2o <- subset(resod2o, abs(log2FoldChange) >=1)
    resod3o <- subset(resod3o, abs(log2FoldChange) >=1)
    
    resod12test <- subset(resod12test, abs(log2FoldChange) >=1)
  #Merging DESEQ results + normalized count data (mRNA)
    #1. merge res (DESEQ results) + dds (normalized count data from DESeqDataSet)
      ##converts res-->df 
      ##converts dds --> df, ensures that we are working w/ normalized counts (normalized=TRUE))
      ##we want to merge dfs according to gene name (by = row.names)
    DEdatareproductive <- merge(as.data.frame(resreproductive), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod12 <- merge(as.data.frame(resod12), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod13 <- merge(as.data.frame(resod13), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod1o <- merge(as.data.frame(resod1o), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod23 <- merge(as.data.frame(resod23), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod2o <- merge(as.data.frame(resod2o), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    DEdataod3o <- merge(as.data.frame(resod3o), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    
    DEdataod12test <- merge(as.data.frame(resod12test), as.data.frame(counts(ddsreproductive, normalized=TRUE)), by="row.names", sort=FALSE)
    #2. Change 1st Column Name = "Gene" 
    names(DEdatareproductive)[1] <- "Gene"
    names(DEdataod12)[1] <- "Gene"
    names(DEdataod13)[1] <- "Gene"
    names(DEdataod1o)[1] <- "Gene"
    names(DEdataod23)[1] <- "Gene"
    names(DEdataod2o)[1] <- "Gene"
    names(DEdataod3o)[1] <- "Gene"
    
    names(DEdataod12test)[1] <- "Gene"
    #head(DEdataod12)
    #3. Create subdfs of each pairwise comparison 
      #a. Get column indices with names containing "od1" or "od2"
    #selected_columns_od12 <- grepl("od1|od2|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
    #selected_columns_od13 <- grepl("od1|od3|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
    #selected_columns_od1o <- grepl("od1|ovary|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
    #selected_columns_od23 <- grepl("od2|od3|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
    #selected_columns_od2o <- grepl("od2|ovary|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
    #selected_columns_od3o <- grepl("od3|ovary|Gene|baseMean|log2FoldChange|lfcSE|stat|pvalue|padj", names(DEdatareproductive))
     #b. Subset the data frame based on selected columns 
    #DEdataod12 <- DEdatareproductive[, selected_columns_od12]
    #DEdataod13 <- DEdatareproductive[, selected_columns_od13]
    #DEdataod1o <- DEdatareproductive[, selected_columns_od1o]
    #DEdataod23 <- DEdatareproductive[, selected_columns_od23]
    #DEdataod2o <- DEdatareproductive[, selected_columns_od2o]
    #DEdataod3o <- DEdatareproductive[, selected_columns_od3o]
    #4. Write results
    write.csv(DEdatareproductive, file="[global]DESEQ-results3.csv")
    write.csv(DEdataod12, file="[od12]DESEQ-results3.csv")
    write.csv(DEdataod13, file="[od13]DESEQ-results3.csv")
    write.csv(DEdataod1o, file="[od1o]DESEQ-results3.csv")
    write.csv(DEdataod23, file="[od23]DESEQ-results3.csv")
    write.csv(DEdataod2o, file="[od2o]DESEQ-results3.csv")
    write.csv(DEdataod3o, file="[od3o]DESEQ-results3.csv")
  
    
    
    

    
  