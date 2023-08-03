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

#Checking to make sure that the # of rows dropped after filtering
nrow(od1_od2)
nrow(od1_od2sig)
nrow(od1_od3)
nrow(od1_od3sig)
nrow(od1_ovary)
nrow(od1_ovarysig)
nrow(od2_od3)
nrow(od2_od3sig)
nrow(od2_ovary)
nrow(od2_ovarysig)
nrow(od3_ovary)
nrow(od3_ovarysig)

##finding the DMR count for each tissue comparison 
nrow(od1_od2sig)
nrow(od1_od3sig)
nrow(od1_ovarysig)
nrow(od2_od3sig)
nrow(od2_ovarysig)
nrow(od3_ovarysig)

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
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052532.1'] <-1
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052533.1'] <-2
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052534.1'] <-3
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052535.1'] <-4
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052536.1'] <-5
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052537.1'] <-6
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052538.1'] <-7
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052539.1'] <-8
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052540.1'] <-9
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052541.1'] <-10
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052542.1'] <-11
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052543.1'] <-12
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052544.1'] <-13
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052545.1'] <-14
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052546.1'] <-15
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052547.1'] <-16
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052548.1'] <-17
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052549.1'] <-18
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052550.1'] <-19
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052551.1'] <-20
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052552.1'] <-21
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052553.1'] <-22
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052554.1'] <-23
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052555.1'] <-24
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052556.1'] <-25
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052557.1'] <-26
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052558.1'] <-27
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052559.1'] <-28
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052560.1'] <-29
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052561.1'] <-30
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052562.1'] <-31
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052563.1'] <-32
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052564.1'] <-33
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052565.1'] <-34
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052566.1'] <-35
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052567.1'] <-36
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052568.1'] <-37
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052569.1'] <-38
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052570.1'] <-39
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052571.1'] <-'W'
od1_od2sig$chr [ od1_od2sig$chr == 'NC_052572.1'] <-'Z'
od1_od2sig$chr [ od1_od2sig$chr == 'NC_053523.1'] <-'MT'
#OD13
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052532.1'] <-1
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052533.1'] <-2
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052534.1'] <-3
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052535.1'] <-4
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052536.1'] <-5
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052537.1'] <-6
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052538.1'] <-7
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052539.1'] <-8
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052540.1'] <-9
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052541.1'] <-10
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052542.1'] <-11
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052543.1'] <-12
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052544.1'] <-13
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052545.1'] <-14
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052546.1'] <-15
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052547.1'] <-16
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052548.1'] <-17
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052549.1'] <-18
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052550.1'] <-19
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052551.1'] <-20
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052552.1'] <-21
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052553.1'] <-22
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052554.1'] <-23
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052555.1'] <-24
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052556.1'] <-25
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052557.1'] <-26
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052558.1'] <-27
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052559.1'] <-28
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052560.1'] <-29
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052561.1'] <-30
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052562.1'] <-31
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052563.1'] <-32
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052564.1'] <-33
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052565.1'] <-34
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052566.1'] <-35
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052567.1'] <-36
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052568.1'] <-37
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052569.1'] <-38
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052570.1'] <-39
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052571.1'] <-'W'
od1_od3sig$chr [ od1_od3sig$chr == 'NC_052572.1'] <-'Z'
od1_od3sig$chr [ od1_od3sig$chr == 'NC_053523.1'] <-'MT'
#OD1o
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052532.1'] <-1
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052533.1'] <-2
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052534.1'] <-3
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052535.1'] <-4
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052536.1'] <-5
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052537.1'] <-6
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052538.1'] <-7
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052539.1'] <-8
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052540.1'] <-9
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052541.1'] <-10
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052542.1'] <-11
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052543.1'] <-12
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052544.1'] <-13
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052545.1'] <-14
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052546.1'] <-15
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052547.1'] <-16
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052548.1'] <-17
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052549.1'] <-18
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052550.1'] <-19
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052551.1'] <-20
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052552.1'] <-21
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052553.1'] <-22
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052554.1'] <-23
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052555.1'] <-24
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052556.1'] <-25
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052557.1'] <-26
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052558.1'] <-27
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052559.1'] <-28
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052560.1'] <-29
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052561.1'] <-30
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052562.1'] <-31
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052563.1'] <-32
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052564.1'] <-33
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052565.1'] <-34
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052566.1'] <-35
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052567.1'] <-36
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052568.1'] <-37
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052569.1'] <-38
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052570.1'] <-39
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052571.1'] <-'W'
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_052572.1'] <-'Z'
od1_ovarysig$chr [ od1_ovarysig$chr == 'NC_053523.1'] <-'MT'
#OD23
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052532.1'] <-1
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052533.1'] <-2
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052534.1'] <-3
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052535.1'] <-4
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052536.1'] <-5
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052537.1'] <-6
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052538.1'] <-7
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052539.1'] <-8
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052540.1'] <-9
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052541.1'] <-10
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052542.1'] <-11
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052543.1'] <-12
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052544.1'] <-13
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052545.1'] <-14
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052546.1'] <-15
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052547.1'] <-16
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052548.1'] <-17
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052549.1'] <-18
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052550.1'] <-19
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052551.1'] <-20
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052552.1'] <-21
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052553.1'] <-22
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052554.1'] <-23
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052555.1'] <-24
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052556.1'] <-25
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052557.1'] <-26
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052558.1'] <-27
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052559.1'] <-28
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052560.1'] <-29
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052561.1'] <-30
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052562.1'] <-31
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052563.1'] <-32
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052564.1'] <-33
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052565.1'] <-34
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052566.1'] <-35
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052567.1'] <-36
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052568.1'] <-37
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052569.1'] <-38
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052570.1'] <-39
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052571.1'] <-'W'
od2_od3sig$chr [ od2_od3sig$chr == 'NC_052572.1'] <-'Z'
od2_od3sig$chr [ od2_od3sig$chr == 'NC_053523.1'] <-'MT'
#OD2o
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052532.1'] <-1
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052533.1'] <-2
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052534.1'] <-3
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052535.1'] <-4
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052536.1'] <-5
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052537.1'] <-6
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052538.1'] <-7
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052539.1'] <-8
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052540.1'] <-9
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052541.1'] <-10
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052542.1'] <-11
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052543.1'] <-12
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052544.1'] <-13
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052545.1'] <-14
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052546.1'] <-15
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052547.1'] <-16
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052548.1'] <-17
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052549.1'] <-18
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052550.1'] <-19
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052551.1'] <-20
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052552.1'] <-21
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052553.1'] <-22
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052554.1'] <-23
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052555.1'] <-24
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052556.1'] <-25
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052557.1'] <-26
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052558.1'] <-27
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052559.1'] <-28
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052560.1'] <-29
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052561.1'] <-30
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052562.1'] <-31
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052563.1'] <-32
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052564.1'] <-33
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052565.1'] <-34
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052566.1'] <-35
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052567.1'] <-36
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052568.1'] <-37
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052569.1'] <-38
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052570.1'] <-39
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052571.1'] <-'W'
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_052572.1'] <-'Z'
od2_ovarysig$chr [ od2_ovarysig$chr == 'NC_053523.1'] <-'MT'
#OD3o
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052532.1'] <-1
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052533.1'] <-2
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052534.1'] <-3
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052535.1'] <-4
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052536.1'] <-5
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052537.1'] <-6
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052538.1'] <-7
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052539.1'] <-8
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052540.1'] <-9
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052541.1'] <-10
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052542.1'] <-11
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052543.1'] <-12
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052544.1'] <-13
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052545.1'] <-14
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052546.1'] <-15
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052547.1'] <-16
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052548.1'] <-17
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052549.1'] <-18
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052550.1'] <-19
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052551.1'] <-20
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052552.1'] <-21
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052553.1'] <-22
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052554.1'] <-23
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052555.1'] <-24
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052556.1'] <-25
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052557.1'] <-26
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052558.1'] <-27
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052559.1'] <-28
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052560.1'] <-29
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052561.1'] <-30
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052562.1'] <-31
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052563.1'] <-32
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052564.1'] <-33
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052565.1'] <-34
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052566.1'] <-35
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052567.1'] <-36
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052568.1'] <-37
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052569.1'] <-38
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052570.1'] <-39
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052571.1'] <-'W'
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_052572.1'] <-'Z'
od3_ovarysig$chr [ od3_ovarysig$chr == 'NC_053523.1'] <-'MT'

#make new df that only isolates our columns of interest (CHR, start, end)
keeps <- c("chr","start","end")
od12= od1_od2sig[keeps]

od13= od1_od3sig[keeps]
od1o= od1_ovarysig[keeps]
od23= od2_od3sig[keeps]
od2o= od2_ovarysig[keeps]
od3o= od3_ovarysig[keeps]

#turning new df into a list (that is compatible w/ biomaRt)
od12list=unite(od12, col="chr-start-end",sep=":")
head(od12list)

od13list=unite(od13, col="chr-start-end",sep=":")
od1olist=unite(od1o, col="chr-start-end",sep=":")
od23list=unite(od23, col="chr-start-end",sep=":")
od2olist=unite(od2o, col="chr-start-end",sep=":")
od3olist=unite(od3o, col="chr-start-end",sep=":")

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
od12biomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od12list, mart = ensembl)

od13biomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od13list, mart = ensembl)

od1obiomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od1olist, mart = ensembl)

od23biomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od23list, mart = ensembl)

od2obiomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od2olist, mart = ensembl)

od3obiomart <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name"
), filters = c("chromosomal_region"), values = od3olist, mart = ensembl)


##Creating output table w/ gene names 
write.table(od12biomart, file= "od12biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(od13biomart, file= "od13biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(od1obiomart, file= "od1obiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(od23biomart, file= "od23biomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(od2obiomart, file= "od2obiomart.txt", row.names = F, sep = "\t",
            quote = F)
write.table(od3obiomart, file= "od3obiomart.txt", row.names = F, sep = "\t",
            quote = F)


