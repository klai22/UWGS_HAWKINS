#Kenneth Lai
# setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/bubbleplot")

#load packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("wesanderson")

#import data 
od12 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od12]enrichment_all.csv")
od13 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od13]enrichment_all.csv")
od1o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od1o]enrichment_all.csv")
od23 = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od23]enrichment_all.csv")
od2o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od2o]enrichment_all.csv")
od3o = read.csv("/Users/kennethlai/desktop/HAWKINS/FAANg/RNAseq/EA/2/[od3o]enrichment_all.csv")

#isolate columns of interest 
keeps <- c("Pathway","Fold.Enrichment","nGenes")
od12pre= od12[keeps]
od13pre= od13[keeps]
od1opre= od1o[keeps]
od23pre= od23[keeps]
od2opre= od2o[keeps]
od3opre= od3o[keeps]

#tagging each column so that I know what column came from what sample in the master file 
od12tag <- od12pre %>% rename(Fold.Enrichment.od12 = Fold.Enrichment, nGenesod12 = nGenes)
od13tag <- od13pre %>% rename(Fold.Enrichment.od13 = Fold.Enrichment, nGenesod13 = nGenes)
od1otag <- od1opre %>% rename(Fold.Enrichment.od1o = Fold.Enrichment, nGenesod1o = nGenes)
od23tag <- od23pre %>% rename(Fold.Enrichment.od23 = Fold.Enrichment, nGenesod23 = nGenes)
od2otag <- od2opre %>% rename(Fold.Enrichment.od2o = Fold.Enrichment, nGenesod2o = nGenes)
od3otag <- od3opre %>% rename(Fold.Enrichment.od3o = Fold.Enrichment, nGenesod3o = nGenes)

#Merging data into a single table 
masterenrichment = merge(od12tag,od13tag, by="Pathway", all=TRUE)
masterenrichment = merge(masterenrichment,od1otag, by="Pathway", all=TRUE)
masterenrichment = merge(masterenrichment,od23tag, by="Pathway", all=TRUE)
masterenrichment = merge(masterenrichment,od2otag, by="Pathway", all=TRUE)
masterenrichment = merge(masterenrichment,od3otag, by="Pathway", all=TRUE)

# Seperate FDR & ngenes into seperate tables--> Replace any missing values (NA) with 0 for ngenes table , keep NAs for FDR 
keepsnGenes <- c("Pathway", "nGenesod12","nGenesod13","nGenesod1o", "nGenesod23","nGenesod2o","nGenesod3o")
keepsFDR <- c("Pathway", "Fold.Enrichment.od12","Fold.Enrichment.od13","Fold.Enrichment.od1o", "Fold.Enrichment.od23","Fold.Enrichment.od2o","Fold.Enrichment.od3o")

masterngenes= masterenrichment[keepsnGenes]
masterFDR= masterenrichment[keepsFDR]

masterngenes[is.na(masterngenes)] <- 0


#Reorganizing to be bubble plot friendly 
masterFDR2 <- masterFDR %>% rename(od12 = Fold.Enrichment.od12,od13 = Fold.Enrichment.od13,od1o = Fold.Enrichment.od1o,od23 = Fold.Enrichment.od23,od2o = Fold.Enrichment.od2o,od3o = Fold.Enrichment.od3o)
masterFDR3 = masterFDR2 %>% 
  pivot_longer(cols=c(od12,od13,od1o,od23,od2o,od3o),names_to="tissuescompared",values_to="Fold.Enrichment")

masterngenes2 <- masterngenes %>% rename(od12 = nGenesod12,od13 = nGenesod13,od1o = nGenesod1o,od23 = nGenesod23,od2o = nGenesod2o,od3o = nGenesod3o)
masterngenes3 = masterngenes2 %>% 
  pivot_longer(cols=c(od12,od13,od1o,od23,od2o,od3o),names_to="tissuescompared",values_to="Transcripts Detected")

#combining these 2 tables 
#making a unique ID ("ID")to match b/w data 
masterFDR4 <- masterFDR3 %>% mutate(ID = paste(Pathway, tissuescompared, sep = "_"))
masterngenes4 <- masterngenes3 %>% mutate(ID = paste(Pathway, tissuescompared, sep = "_"))
#merging dfs by unique ID 
masterbubbledata<- merge(masterFDR4, masterngenes4, by = "ID")
#eliminating duplicated columns 
masterbubbledata2 = masterbubbledata %>% select(-one_of(c("ID","Pathway.y","tissuescompared.y")))


#BUBBLE PLOT (size = ngenes, color = FDR)
EAbubble=ggplot(masterbubbledata2)+geom_point(aes(x=tissuescompared.x,Pathway.x,size=`Transcripts Detected`,color=Fold.Enrichment)) +
      labs(x="Tissue Compared", y="KEGG Pathway",color="Fold Enrichment (FDR)",size="Number of Detected Transcripts") +
      theme_linedraw() +
      theme(axis.text.x=element_text(angle=25, vjust=0.88,hjust=0.9),
            plot.margin=margin(0.8,0.8,0.8,0.8,"cm"))+
      scale_color_gradient(low = "blue", high = "red")
    

#printing the graph 
pdf("[1]EAbubble.pdf", width=10,height=8)

print(EAbubble,
      labels=c("Enrichment Analysis"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()
