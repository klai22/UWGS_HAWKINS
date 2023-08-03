##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library (biomaRt)
library(tidyr)

#import data
od1_od2=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od1-od2_filter_DMR_mod.txt",header=TRUE)
od1_od3=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od1-od3_filter_DMR_mod.txt",header=TRUE)
od1_ovary=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od1-ovary_filter_DMR_mod.txt",header=TRUE)
od2_od3=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od2-od3_filter_DMR_mod.txt",header=TRUE)
od2_ovary=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od2-ovary_filter_DMR_mod.txt",header=TRUE)
od3_ovary=read.table("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/od3-ovary_filter_DMR_mod.txt",header=TRUE)

#Remove rows that do not fit criteria (removing rows that have "Inf" in the t-test column + that have threshold values <0.0625 [the new])
od1_od2sig = subset(od1_od2, t != 'Inf' & Thershold >= 0.0625)
od1_od3sig = subset(od1_od3, t != 'Inf' & Thershold >= 0.0625)
od1_ovarysig=subset(od1_ovary, t != 'Inf' & Thershold >= 0.0625)
od2_od3sig=subset(od2_od3, t != 'Inf' & Thershold >= 0.0625)
od2_ovarysig=subset(od2_ovary, t != 'Inf' & Thershold >= 0.0625)
od3_ovarysig=subset(od3_ovary, t != 'Inf' & Thershold >= 0.0625)

#Checking to make sure that the # of rows dropped after filtering
nrow(od1_od2)
nrow(od1_od2sig)

#make new df that only isolates our columns of interest (CHR, start, end)
keeps <- c("chr","start","end")
od12pre= od1_od2sig[keeps]
od13pre= od1_od3sig[keeps]
od1opre= od1_ovarysig[keeps]
od23pre= od2_od3sig[keeps]
od2opre= od2_ovarysig[keeps]
od3opre= od3_ovarysig[keeps]

#turning our dfs --> tab-delimited bed files 
od12 <- unite(od12pre, col = "chr_start_end", sep = "\t")
od13 <- unite(od13pre, col = "chr_start_end", sep = "\t")
od1o <- unite(od1opre, col = "chr_start_end", sep = "\t")
od23 <- unite(od23pre, col = "chr_start_end", sep = "\t")
od2o <- unite(od2opre, col = "chr_start_end", sep = "\t")
od3o <- unite(od3opre, col = "chr_start_end", sep = "\t")

head(od12)
head(od13)
head(od1o)
head(od23)
head(od2o)
head(od3o)


# write the data frame to a tab-separated text file
write.table(od12, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od12homerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
write.table(od13, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od13homerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
write.table(od1o, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od1ohomerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
write.table(od23, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od23homerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
write.table(od2o, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od2ohomerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
write.table(od3o, file = "/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer/od3ohomerupdated.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)



