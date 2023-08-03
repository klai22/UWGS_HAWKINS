##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED PARAMETERS/Enrichment Analysis")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library (biomaRt)
library(tidyr)

#import data
od1_od2=read.table("od1-od2_filter_DMR_mod.txt",header=TRUE, row.names = 1)
od1_od3=read.table("od1-od3_filter_DMR_mod.txt",header=TRUE, row.names = 1)
od1_ovary=read.table("od1-ovary_filter_DMR_mod.txt",header=TRUE, row.names = 1)
od2_od3=read.table("od2-od3_filter_DMR_mod.txt",header=TRUE, row.names = 1)
od2_ovary=read.table("od2-ovary_filter_DMR_mod.txt",header=TRUE, row.names = 1)
od3_ovary=read.table("od3-ovary_filter_DMR_mod.txt",header=TRUE, row.names = 1)

#Remove rows that do not fit criteria (removing rows that have "Inf" in the t-test column + that have threshold values <0.0625 [the new])
od1_od2sig = subset(od1_od2, t != 'Inf' & Thershold >= 0.0625)
od1_od3sig = subset(od1_od3, t != 'Inf' & Thershold >= 0.0625)
od1_ovarysig=subset(od1_ovary, t != 'Inf' & Thershold >= 0.0625)
od2_od3sig=subset(od2_od3, t != 'Inf' & Thershold >= 0.0625)
od2_ovarysig=subset(od2_ovary, t != 'Inf' & Thershold >= 0.0625)
od3_ovarysig=subset(od3_ovary, t != 'Inf' & Thershold >= 0.0625)

#Isolating Hypo vs. Hyper DMRs 
hyperod1_od2 = od1_od2sig[od1_od2sig$t >= 0,]
hypood1_od2 = od1_od2sig[od1_od2sig$t < 0,]
hyperod1_od3 = od1_od3sig[od1_od3sig$t >= 0,]
hypood1_od3 = od1_od3sig[od1_od3sig$t < 0,]
hyperod1_ovary = od1_ovarysig[od1_ovarysig$t >= 0,]
hypood1_ovary = od1_ovarysig[od1_ovarysig$t < 0,]
hyperod2_od3 = od2_od3sig[od2_od3sig$t >= 0,]
hypood2_od3 = od2_od3sig[od2_od3sig$t < 0,]
hyperod2_ovary = od2_ovarysig[od2_ovarysig$t >= 0,]
hypood2_ovary = od2_ovarysig[od2_ovarysig$t < 0,]
hyperod3_ovary = od3_ovarysig[od3_ovarysig$t >= 0,]
hypood3_ovary = od3_ovarysig[od3_ovarysig$t < 0,]


#Checking to make sure that the # of rows dropped after filtering
nrow(od1_od2)
nrow(od1_od2sig)
nrow(hyperod1_od2)
nrow(hypood1_od2)
nrow(od1_od3)
nrow(od1_od3sig)
nrow(hyperod1_od3)
nrow(hypood1_od3)
##skipped the others to save time 

##BIOMART 
#Install 
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")

#fixing install issues
####install.packages("lifecycle")

###Prepping Data for BiomaRt 
##making my dfs compatible w/ biomaRt (replacing chromosome ref seq.s w/ chromosome names : https://www.ncbi.nlm.nih.gov/genome/?term=chicken)
#OD12
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052532.1'] <-1
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052532.1'] <-1
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052533.1'] <-2
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052533.1'] <-2
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052534.1'] <-3
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052534.1'] <-3
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052535.1'] <-4
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052535.1'] <-4
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052536.1'] <-5
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052536.1'] <-5
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052537.1'] <-6
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052537.1'] <-6
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052538.1'] <-7
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052538.1'] <-7
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052539.1'] <-8
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052539.1'] <-8
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052540.1'] <-9
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052540.1'] <-9
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052541.1'] <-10
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052541.1'] <-10
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052542.1'] <-11
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052542.1'] <-11
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052543.1'] <-12
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052543.1'] <-12
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052544.1'] <-13
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052544.1'] <-13
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052545.1'] <-14
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052545.1'] <-14
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052546.1'] <-15
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052546.1'] <-15
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052547.1'] <-16
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052547.1'] <-16
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052548.1'] <-17
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052548.1'] <-17
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052549.1'] <-18
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052549.1'] <-18
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052550.1'] <-19
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052550.1'] <-19
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052551.1'] <-20
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052551.1'] <-20
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052552.1'] <-21
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052552.1'] <-21
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052553.1'] <-22
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052553.1'] <-22
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052554.1'] <-23
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052554.1'] <-23
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052555.1'] <-24
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052555.1'] <-24
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052556.1'] <-25
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052556.1'] <-25
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052557.1'] <-26
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052557.1'] <-26
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052558.1'] <-27
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052558.1'] <-27
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052559.1'] <-28
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052559.1'] <-28
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052560.1'] <-29
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052560.1'] <-29
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052561.1'] <-30
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052561.1'] <-30
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052562.1'] <-31
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052562.1'] <-31
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052563.1'] <-32
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052563.1'] <-32
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052564.1'] <-33
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052564.1'] <-33
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052565.1'] <-34
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052565.1'] <-34
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052566.1'] <-35
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052566.1'] <-35
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052567.1'] <-36
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052567.1'] <-36
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052568.1'] <-37
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052568.1'] <-37
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052569.1'] <-38
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052569.1'] <-38
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052570.1'] <-39
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052570.1'] <-39
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052571.1'] <-'W'
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052571.1'] <-'W'
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_052572.1'] <-'Z'
hypood1_od2$chr [ hypood1_od2$chr == 'NC_052572.1'] <-'Z'
hyperod1_od2$chr [ hyperod1_od2$chr == 'NC_053523.1'] <-'MT'
hypood1_od2$chr [ hypood1_od2$chr == 'NC_053523.1'] <-'MT'
#OD13
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052532.1'] <-1
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052532.1'] <-1
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052533.1'] <-2
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052533.1'] <-2
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052534.1'] <-3
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052534.1'] <-3
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052535.1'] <-4
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052535.1'] <-4
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052536.1'] <-5
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052536.1'] <-5
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052537.1'] <-6
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052537.1'] <-6
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052538.1'] <-7
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052538.1'] <-7
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052539.1'] <-8
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052539.1'] <-8
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052540.1'] <-9
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052540.1'] <-9
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052541.1'] <-10
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052541.1'] <-10
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052542.1'] <-11
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052542.1'] <-11
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052543.1'] <-12
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052543.1'] <-12
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052544.1'] <-13
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052544.1'] <-13
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052545.1'] <-14
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052545.1'] <-14
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052546.1'] <-15
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052546.1'] <-15
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052547.1'] <-16
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052547.1'] <-16
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052548.1'] <-17
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052548.1'] <-17
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052549.1'] <-18
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052549.1'] <-18
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052550.1'] <-19
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052550.1'] <-19
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052551.1'] <-20
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052551.1'] <-20
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052552.1'] <-21
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052552.1'] <-21
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052553.1'] <-22
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052553.1'] <-22
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052554.1'] <-23
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052554.1'] <-23
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052555.1'] <-24
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052555.1'] <-24
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052556.1'] <-25
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052556.1'] <-25
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052557.1'] <-26
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052557.1'] <-26
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052558.1'] <-27
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052558.1'] <-27
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052559.1'] <-28
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052559.1'] <-28
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052560.1'] <-29
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052560.1'] <-29
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052561.1'] <-30
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052561.1'] <-30
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052562.1'] <-31
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052562.1'] <-31
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052563.1'] <-32
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052563.1'] <-32
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052564.1'] <-33
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052564.1'] <-33
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052565.1'] <-34
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052565.1'] <-34
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052566.1'] <-35
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052566.1'] <-35
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052567.1'] <-36
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052567.1'] <-36
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052568.1'] <-37
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052568.1'] <-37
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052569.1'] <-38
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052569.1'] <-38
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052570.1'] <-39
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052570.1'] <-39
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052571.1'] <-'W'
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052571.1'] <-'W'
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_052572.1'] <-'Z'
hypood1_od3$chr [ hypood1_od3$chr == 'NC_052572.1'] <-'Z'
hyperod1_od3$chr [ hyperod1_od3$chr == 'NC_053523.1'] <-'MT'
hypood1_od3$chr [ hypood1_od3$chr == 'NC_053523.1'] <-'MT'
#OD1OVARY 
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052532.1'] <-1
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052532.1'] <-1
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052533.1'] <-2
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052533.1'] <-2
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052534.1'] <-3
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052534.1'] <-3
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052535.1'] <-4
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052535.1'] <-4
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052536.1'] <-5
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052536.1'] <-5
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052537.1'] <-6
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052537.1'] <-6
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052538.1'] <-7
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052538.1'] <-7
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052539.1'] <-8
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052539.1'] <-8
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052540.1'] <-9
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052540.1'] <-9
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052541.1'] <-10
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052541.1'] <-10
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052542.1'] <-11
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052542.1'] <-11
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052543.1'] <-12
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052543.1'] <-12
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052544.1'] <-13
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052544.1'] <-13
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052545.1'] <-14
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052545.1'] <-14
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052546.1'] <-15
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052546.1'] <-15
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052547.1'] <-16
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052547.1'] <-16
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052548.1'] <-17
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052548.1'] <-17
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052549.1'] <-18
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052549.1'] <-18
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052550.1'] <-19
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052550.1'] <-19
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052551.1'] <-20
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052551.1'] <-20
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052552.1'] <-21
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052552.1'] <-21
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052553.1'] <-22
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052553.1'] <-22
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052554.1'] <-23
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052554.1'] <-23
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052555.1'] <-24
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052555.1'] <-24
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052556.1'] <-25
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052556.1'] <-25
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052557.1'] <-26
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052557.1'] <-26
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052558.1'] <-27
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052558.1'] <-27
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052559.1'] <-28
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052559.1'] <-28
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052560.1'] <-29
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052560.1'] <-29
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052561.1'] <-30
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052561.1'] <-30
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052562.1'] <-31
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052562.1'] <-31
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052563.1'] <-32
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052563.1'] <-32
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052564.1'] <-33
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052564.1'] <-33
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052565.1'] <-34
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052565.1'] <-34
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052566.1'] <-35
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052566.1'] <-35
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052567.1'] <-36
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052567.1'] <-36
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052568.1'] <-37
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052568.1'] <-37
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052569.1'] <-38
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052569.1'] <-38
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052570.1'] <-39
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052570.1'] <-39
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052571.1'] <-'W'
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052571.1'] <-'W'
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_052572.1'] <-'Z'
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_052572.1'] <-'Z'
hyperod1_ovary$chr [ hyperod1_ovary$chr == 'NC_053523.1'] <-'MT'
hypood1_ovary$chr [ hypood1_ovary$chr == 'NC_053523.1'] <-'MT'
#OD23
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052532.1'] <-1
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052532.1'] <-1
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052533.1'] <-2
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052533.1'] <-2
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052534.1'] <-3
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052534.1'] <-3
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052535.1'] <-4
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052535.1'] <-4
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052536.1'] <-5
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052536.1'] <-5
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052537.1'] <-6
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052537.1'] <-6
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052538.1'] <-7
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052538.1'] <-7
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052539.1'] <-8
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052539.1'] <-8
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052540.1'] <-9
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052540.1'] <-9
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052541.1'] <-10
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052541.1'] <-10
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052542.1'] <-11
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052542.1'] <-11
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052543.1'] <-12
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052543.1'] <-12
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052544.1'] <-13
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052544.1'] <-13
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052545.1'] <-14
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052545.1'] <-14
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052546.1'] <-15
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052546.1'] <-15
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052547.1'] <-16
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052547.1'] <-16
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052548.1'] <-17
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052548.1'] <-17
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052549.1'] <-18
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052549.1'] <-18
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052550.1'] <-19
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052550.1'] <-19
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052551.1'] <-20
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052551.1'] <-20
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052552.1'] <-21
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052552.1'] <-21
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052553.1'] <-22
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052553.1'] <-22
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052554.1'] <-23
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052554.1'] <-23
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052555.1'] <-24
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052555.1'] <-24
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052556.1'] <-25
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052556.1'] <-25
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052557.1'] <-26
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052557.1'] <-26
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052558.1'] <-27
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052558.1'] <-27
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052559.1'] <-28
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052559.1'] <-28
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052560.1'] <-29
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052560.1'] <-29
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052561.1'] <-30
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052561.1'] <-30
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052562.1'] <-31
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052562.1'] <-31
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052563.1'] <-32
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052563.1'] <-32
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052564.1'] <-33
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052564.1'] <-33
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052565.1'] <-34
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052565.1'] <-34
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052566.1'] <-35
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052566.1'] <-35
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052567.1'] <-36
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052567.1'] <-36
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052568.1'] <-37
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052568.1'] <-37
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052569.1'] <-38
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052569.1'] <-38
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052570.1'] <-39
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052570.1'] <-39
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052571.1'] <-'W'
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052571.1'] <-'W'
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_052572.1'] <-'Z'
hypood2_od3$chr [ hypood2_od3$chr == 'NC_052572.1'] <-'Z'
hyperod2_od3$chr [ hyperod2_od3$chr == 'NC_053523.1'] <-'MT'
hypood2_od3$chr [ hypood2_od3$chr == 'NC_053523.1'] <-'MT'
#OD2OVARY
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052532.1'] <-1
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052532.1'] <-1
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052533.1'] <-2
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052533.1'] <-2
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052534.1'] <-3
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052534.1'] <-3
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052535.1'] <-4
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052535.1'] <-4
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052536.1'] <-5
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052536.1'] <-5
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052537.1'] <-6
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052537.1'] <-6
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052538.1'] <-7
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052538.1'] <-7
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052539.1'] <-8
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052539.1'] <-8
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052540.1'] <-9
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052540.1'] <-9
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052541.1'] <-10
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052541.1'] <-10
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052542.1'] <-11
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052542.1'] <-11
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052543.1'] <-12
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052543.1'] <-12
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052544.1'] <-13
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052544.1'] <-13
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052545.1'] <-14
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052545.1'] <-14
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052546.1'] <-15
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052546.1'] <-15
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052547.1'] <-16
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052547.1'] <-16
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052548.1'] <-17
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052548.1'] <-17
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052549.1'] <-18
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052549.1'] <-18
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052550.1'] <-19
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052550.1'] <-19
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052551.1'] <-20
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052551.1'] <-20
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052552.1'] <-21
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052552.1'] <-21
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052553.1'] <-22
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052553.1'] <-22
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052554.1'] <-23
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052554.1'] <-23
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052555.1'] <-24
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052555.1'] <-24
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052556.1'] <-25
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052556.1'] <-25
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052557.1'] <-26
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052557.1'] <-26
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052558.1'] <-27
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052558.1'] <-27
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052559.1'] <-28
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052559.1'] <-28
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052560.1'] <-29
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052560.1'] <-29
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052561.1'] <-30
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052561.1'] <-30
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052562.1'] <-31
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052562.1'] <-31
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052563.1'] <-32
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052563.1'] <-32
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052564.1'] <-33
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052564.1'] <-33
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052565.1'] <-34
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052565.1'] <-34
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052566.1'] <-35
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052566.1'] <-35
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052567.1'] <-36
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052567.1'] <-36
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052568.1'] <-37
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052568.1'] <-37
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052569.1'] <-38
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052569.1'] <-38
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052570.1'] <-39
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052570.1'] <-39
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052571.1'] <-'W'
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052571.1'] <-'W'
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_052572.1'] <-'Z'
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_052572.1'] <-'Z'
hyperod2_ovary$chr [ hyperod2_ovary$chr == 'NC_053523.1'] <-'MT'
hypood2_ovary$chr [ hypood2_ovary$chr == 'NC_053523.1'] <-'MT'
#OD3OVARY 
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052532.1'] <-1
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052532.1'] <-1
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052533.1'] <-2
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052533.1'] <-2
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052534.1'] <-3
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052534.1'] <-3
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052535.1'] <-4
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052535.1'] <-4
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052536.1'] <-5
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052536.1'] <-5
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052537.1'] <-6
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052537.1'] <-6
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052538.1'] <-7
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052538.1'] <-7
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052539.1'] <-8
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052539.1'] <-8
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052540.1'] <-9
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052540.1'] <-9
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052541.1'] <-10
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052541.1'] <-10
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052542.1'] <-11
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052542.1'] <-11
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052543.1'] <-12
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052543.1'] <-12
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052544.1'] <-13
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052544.1'] <-13
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052545.1'] <-14
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052545.1'] <-14
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052546.1'] <-15
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052546.1'] <-15
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052547.1'] <-16
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052547.1'] <-16
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052548.1'] <-17
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052548.1'] <-17
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052549.1'] <-18
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052549.1'] <-18
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052550.1'] <-19
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052550.1'] <-19
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052551.1'] <-20
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052551.1'] <-20
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052552.1'] <-21
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052552.1'] <-21
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052553.1'] <-22
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052553.1'] <-22
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052554.1'] <-23
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052554.1'] <-23
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052555.1'] <-24
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052555.1'] <-24
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052556.1'] <-25
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052556.1'] <-25
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052557.1'] <-26
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052557.1'] <-26
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052558.1'] <-27
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052558.1'] <-27
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052559.1'] <-28
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052559.1'] <-28
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052560.1'] <-29
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052560.1'] <-29
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052561.1'] <-30
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052561.1'] <-30
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052562.1'] <-31
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052562.1'] <-31
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052563.1'] <-32
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052563.1'] <-32
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052564.1'] <-33
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052564.1'] <-33
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052565.1'] <-34
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052565.1'] <-34
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052566.1'] <-35
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052566.1'] <-35
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052567.1'] <-36
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052567.1'] <-36
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052568.1'] <-37
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052568.1'] <-37
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052569.1'] <-38
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052569.1'] <-38
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052570.1'] <-39
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052570.1'] <-39
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052571.1'] <-'W'
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052571.1'] <-'W'
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_052572.1'] <-'Z'
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_052572.1'] <-'Z'
hyperod3_ovary$chr [ hyperod3_ovary$chr == 'NC_053523.1'] <-'MT'
hypood3_ovary$chr [ hypood3_ovary$chr == 'NC_053523.1'] <-'MT'


#make new df that only isolates our columns of interest (CHR, start, end)
keeps <- c("chr","start","end")
hyperod12= hyperod1_od2[keeps]
hypood12=hypood1_od2[keeps]
hyperod13= hyperod1_od3[keeps]
hypood13=hypood1_od3[keeps]
hyperod1ovary= hyperod1_ovary[keeps]
hypood1ovary=hypood1_ovary[keeps]
hyperod23= hyperod2_od3[keeps]
hypood23=hypood2_od3[keeps]
hyperod2ovary= hyperod2_ovary[keeps]
hypood2ovary=hypood2_ovary[keeps]
hyperod3ovary= hyperod3_ovary[keeps]
hypood3ovary=hypood3_ovary[keeps]


#turning new df into a list (that is compatible w/ biomaRt)
hyperod12list=unite(hyperod12, col="chr-start-end",sep=":")
hypood12list=unite(hypood12, col="chr-start-end",sep=":")
hyperod13list=unite(hyperod13, col="chr-start-end",sep=":")
hypood13list=unite(hypood13, col="chr-start-end",sep=":")
hyperod1ovarylist=unite(hyperod1ovary, col="chr-start-end",sep=":")
hypood1ovarylist=unite(hypood1ovary, col="chr-start-end",sep=":")
hyperod23list=unite(hyperod23, col="chr-start-end",sep=":")
hypood23list=unite(hypood23, col="chr-start-end",sep=":")
hyperod2ovarylist=unite(hyperod2ovary, col="chr-start-end",sep=":")
hypood2ovarylist=unite(hypood2ovary, col="chr-start-end",sep=":")
hyperod3ovarylist=unite(hyperod3ovary, col="chr-start-end",sep=":")
hypood3ovarylist=unite(hypood3ovary, col="chr-start-end",sep=":")


head(hyperod12list)
head(hypood12list)


###Running BiomaRT(https://bio340files.s3.amazonaws.com/SearchBiomart.html)
#list database 
myMarts <- listMarts()
#select database
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
#list datasets for given database
myDatasets<-listDatasets(ensembl)
#select dataset 
ensembl=useMart("ENSEMBL_MART_ENSEMBL",
                dataset="ggallus_gene_ensembl")
#list filters (telling biomaRt the type of data your input file that aligns/shares w/ the database so that biomaRt can properly match corresponding rows)
filters = listFilters (ensembl)
head(filters)
#List Attributes (attributes = the columns / data you're gonna download from database w/ corresponding data)
atrributes = listAttributes (ensembl)
head(atrributes)
###Running BiomaRt while specifying variables for (1) attributes + (2) filters 
hyperod12biomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod12list, mart = ensembl)

hypood12biomart <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"
), filters = c("chromosomal_region"), values = hypood12list, mart = ensembl)

hyperod13biomart <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod13list, mart = ensembl)

hypood13biomart <- getBM(attributes=c("ensembl_gene_id",
                                      "external_gene_name"
), filters = c("chromosomal_region"), values = hypood13list, mart = ensembl)

hyperod1ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod1ovarylist, mart = ensembl)

hypood1ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                      "external_gene_name"
), filters = c("chromosomal_region"), values = hypood1ovarylist, mart = ensembl)

hyperod23biomart <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod23list, mart = ensembl)

hypood23biomart <- getBM(attributes=c("ensembl_gene_id",
                                      "external_gene_name"
), filters = c("chromosomal_region"), values = hypood23list, mart = ensembl)

hyperod2ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                           "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod2ovarylist, mart = ensembl)

hypood2ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                          "external_gene_name"
), filters = c("chromosomal_region"), values = hypood2ovarylist, mart = ensembl)

hyperod3ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                           "external_gene_name"
), filters = c("chromosomal_region"), values = hyperod3ovarylist, mart = ensembl)

hypood3ovarybiomart <- getBM(attributes=c("ensembl_gene_id",
                                          "external_gene_name"
), filters = c("chromosomal_region"), values = hypood3ovarylist, mart = ensembl)

##Creating output table w/ gene names 
write.table(hyperod12biomart, file= "hyperod12biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood12biomart, file= "hypood12biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hyperod13biomart, file= "hyperod13biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood13biomart, file= "hypood13biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hyperod1ovarybiomart, file= "hyperod1ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood1ovarybiomart, file= "hypood1ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hyperod23biomart, file= "hyperod23biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood23biomart, file= "hypood23biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hyperod2ovarybiomart, file= "hyperod2ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood2ovarybiomart, file= "hypood2ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hyperod3ovarybiomart, file= "hyperod3ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(hypood3ovarybiomart, file= "hypood3ovarybiomart.txt", row.names = F, sep = "\t",
            quote = F)






