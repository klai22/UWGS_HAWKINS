# setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Enrichment Analysis/ShinyGO Results")

#load packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")
library(lubridate)

#import data 
od12 = read.csv("hypood12enrichment.csv")
od13 = read.csv("hypood13enrichment.csv")
od1o = read.csv("hypood1oenrichment.csv")
od23 = read.csv("hypood23enrichment.csv")
od2o = read.csv("hypood2oenrichment.csv")
od3o = read.csv("hypood3oenrichment.csv")

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


###FDR HEATMAP 
#Reorganizing data to be heatmap friendly 
#renaming columns 
masterFDR2 <- masterFDR %>% rename(od12 = Fold.Enrichment.od12,od13 = Fold.Enrichment.od13,od1o = Fold.Enrichment.od1o,od23 = Fold.Enrichment.od23,od2o = Fold.Enrichment.od2o,od3o = Fold.Enrichment.od3o)
#moving fold ernichment data arround to make 2 columns out of it so that it can be used in heatmap 
masterFDR3 = masterFDR2 %>% 
  pivot_longer(cols=c(od12,od13,od1o,od23,od2o,od3o),names_to="tissuescompared",values_to="FDR.Enrichment")

#finding the median value of FDRenrichment so I can apply it to the color correclty 
median(masterFDR3$FDR.Enrichment, na.rm = TRUE)
# the median = 1.465989

#finding the average as a replacement for the median value 
mean(range(masterFDR3$FDR.Enrichment,na.rm = TRUE))
# average = 2.01115

#ggplot for the heatmap 
FDRheatmap=ggplot(masterFDR3)+ geom_tile(aes(x=tissuescompared,y=Pathway,fill=FDR.Enrichment))+  scale_colour_gradient2(
                                                                               low ="lightblue",
                                                                               mid = "black",
                                                                               high = "pink",
                                                                               midpoint = 2.01115,
                                                                               space = "Lab",
                                                                               na.value = "grey50",
                                                                               guide = "colourbar",
                                                                               aesthetics = "fill")+
  labs(x="Tissues Compared", y="KEGG Pathway",fill="FDR Enrichment Value") + 
  scale_x_discrete(labels=c("od12"="od1 vs. od2",
                            "od13"="od1 vs. od3",
                            "od1o"="od1 vs. ovary",
                            "od23"="od2 vs. od3",
                            "od2o"="od2 vs. ovary",
                            "od3o"="od3 vs. ovary"))+
  theme(axis.text.x=element_text(angle=20,size=9, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=7),
        axis.title.x=element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.margin=margin(0.8,0.8,0.8,0.8,"cm"),legend.position="right")+ ggtitle("Hypomethylated DMR KEGG Pathways (FDR Enrichment-Based)")

#printing the graph 
pdf("HypoFDR.pdf", width=10,height=8)

print(FDRheatmap,
      labels=c("Hypomethylated DMR KEGG Pathways (FDR Enrichment-Based)"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()

###NGENE HEATMAP 
#Reorganizing data to be heatmap friendly 
#renaming columns 
masterngenes2 <- masterngenes %>% rename(od12 = nGenesod12,od13 = nGenesod13,od1o = nGenesod1o,od23 = nGenesod23,od2o = nGenesod2o,od3o = nGenesod3o)
#moving fold ernichment data arround to make 2 columns out of it so that it can be used in heatmap 
masterngenes3 = masterngenes2 %>% 
  pivot_longer(cols=c(od12,od13,od1o,od23,od2o,od3o),names_to="tissuescompared",values_to="DMRs Detected")

#finding the median value of FDRenrichment so I can apply it to the color correclty 
median(masterngenes3$`DMRs Detected`)
# the median = 25

#finding the average as a replacement for the median value 
mean(range(masterngenes3$`DMRs Detected`))
# average = 87.5

##!!! UNLIKE FDR ENRICHMENT, nGENES is a DISCRETE VARIABLE (not continous)so we have to make it a factor variable to be plotted on a heatmap, and also edit the heatmap code a little 
# Convert DMRs Detected to a factor variable
#masterngenes3$"DMRs Detected"<- factor(masterngenes3$"DMRs Detected")

#I tried plotting DMRs detected as a discrete variable and it fave me all blocks of the same color, so for now I will treat it as a continous variable...

#ggplot for the heatmap 
nDMRsheatmap=ggplot(masterngenes3)+ geom_tile(aes(x=tissuescompared,y=Pathway,fill=`DMRs Detected`))+ scale_colour_gradient2(
  low ="lightblue",
  mid = "black",
  high = "pink",
  midpoint = 87.5,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill")+
  labs(x="tissuescompared", y="KEGG Pathway",fill="DMRs Detected") + 
  scale_x_discrete(labels=c("od12"="Magnum vs. Isthmus(od12)",
                            "od13"="Magnum vs. Shell Gland(od13)",
                            "od1o"="Magnum vs. Ovary(od1o)",
                            "od23"="Isthmus vs. Shell Gland (od23)",
                            "od2o"="Isthmus vs. Ovary(od2o)",
                            "od3o"="Shell Gland vs. Ovary(od3o)"))+
  theme(axis.text.x=element_text(angle=30,size=9, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=7),
        axis.title.x=element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.margin=margin(0.8,0.8,0.8,0.8,"cm"),legend.position="right")+ ggtitle("Hypomethylated DMR KEGG Pathways (nGenes-based)")


#printing the graph 
pdf("HyponDMR.pdf", width=10,height=8)

print(nDMRsheatmap,
      labels=c("Hypomethylated DMR KEGG Pathways (nGenes-based)"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()







