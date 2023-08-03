##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/CpGContext")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
require(RColorBrewer)
library(tibble)

##import data
od1=read.csv("od1.mstat.csv",header=TRUE)
od2=read.csv("od2.mstat.csv",header=TRUE)
od3=read.csv("od3.mstat.csv",header=TRUE)
ovary=read.csv("ovary.mstat.csv",header=TRUE)

##CHROMOSOME vs. CG METHYLATION
##Isolated columns + rows of interest
rmrod1=od1[-c(1:10), ]
rmcod1=rmrod1[-c(1,3,5:11)]

rmrod2=od2[-c(1:10), ]
rmcod2=rmrod2[-c(1,3,5:11)]

rmrod3=od3[-c(1:10), ]
rmcod3=rmrod3[-c(1,3,5:11)]

rmrovary=ovary[-c(1:10), ]
rmcovary=rmrovary[-c(1,3,5:11)]

##make the bar graph 
CGmetod1=ggplot(rmcod1,aes(x=context,y=CG)) + 
  geom_col(fill="pink",width=0.7) + 
  xlab("Chromosome")+
  ylab("CG Methylation Lvl") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD1 CG Methylation")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=4, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

CGmetod2=ggplot(rmcod2,aes(x=context,y=CG)) + 
  geom_col(fill="yellow",width=0.7) + 
  xlab("Chromosome")+
  ylab("CG Methylation Lvl") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD2 CG Methylation")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=4, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

CGmetod3=ggplot(rmcod3,aes(x=context,y=CG)) + 
  geom_col(fill="lightgreen",width=0.7) + 
  xlab("Chromosome")+
  ylab("CG Methylation Lvl") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD3 CG Methylation")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=4, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

CGmetovary=ggplot(rmcovary,aes(x=context,y=CG)) + 
  geom_col(fill="skyblue",width=0.7) + 
  xlab("Chromosome")+
  ylab("CG Methylation Lvl") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OVARY CG Methylation")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=4, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

pdf("CGmetChromosomeDistribution.pdf", width=8,height=15)
print(ggarrange(CGmetod1,CGmetod2,CGmetod3,CGmetovary,ncol=1,nrow=4))
dev.off()

##Methylation-Type Comparison 
##Isolated columns + rows of interest
rmrod1A=od1[-c(2:52), ]
rmcod1A=rmrod1A[-c(1:3)]

rmrod2A=od2[-c(2:52), ]
rmcod2A=rmrod2A[-c(1:3)]

rmrod3A=od3[-c(2:52), ]
rmcod3A=rmrod3A[-c(1:3)]

rmrovaryA=ovary[-c(2:52), ]
rmcovaryA=rmrovaryA[-c(1:3)]

##flipping table 
flipod1.a=t(rmcod1A)
flipod1.b=as.data.frame(flipod1.a)
flipod1.c=tibble::rownames_to_column(flipod1.b)
colnames(flipod1.c)=c("nucleotide pair type" , "avg met lvl")

flipod2.a=t(rmcod2A)
flipod2.b=as.data.frame(flipod2.a)
flipod2.c=tibble::rownames_to_column(flipod2.b)
colnames(flipod2.c)=c("nucleotide pair type" , "avg met lvl")

flipod3.a=t(rmcod3A)
flipod3.b=as.data.frame(flipod3.a)
flipod3.c=tibble::rownames_to_column(flipod3.b)
colnames(flipod3.c)=c("nucleotide pair type" , "avg met lvl")

flipovary.a=t(rmcovaryA)
flipovary.b=as.data.frame(flipovary.a)
flipovary.c=tibble::rownames_to_column(flipovary.b)
colnames(flipovary.c)=c("nucleotide pair type" , "avg met lvl")

##make the pie chart..... I GIVE UP, made a bar graph instead....
MetTypeod1=ggplot(flipod1.c,aes(x=flipod1.c$`nucleotide pair type`,y=flipod1.c$`avg met lvl`)) + 
  geom_col(fill="blue",width=0.7) + 
  xlab("Nudelotide Pair Type ")+
  ylab("Average Methylation") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD1 Methylation Distribution")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=6, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

MetTypeod2=ggplot(flipod2.c,aes(x=flipod2.c$`nucleotide pair type`,y=flipod2.c$`avg met lvl`)) + 
  geom_col(fill="blue",width=0.7) + 
  xlab("Nudelotide Pair Type ")+
  ylab("Average Methylation") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD2 Methylation Distribution")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=6, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

MetTypeod3=ggplot(flipod3.c,aes(x=flipod3.c$`nucleotide pair type`,y=flipod3.c$`avg met lvl`)) + 
  geom_col(fill="blue",width=0.7) + 
  xlab("Nudelotide Pair Type ")+
  ylab("Average Methylation") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("OD3 Methylation Distribution")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=6, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))

MetTypeovary=ggplot(flipovary.c,aes(x=flipovary.c$`nucleotide pair type`,y=flipovary.c$`avg met lvl`)) + 
  geom_col(fill="blue",width=0.7) + 
  xlab("Nudelotide Pair Type ")+
  ylab("Average Methylation") +
  ylim(0,1)+
  theme_classic() +
  ggtitle("Ovary Methylation Distribution")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=6, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.margin=margin(1,1,1,1,"cm"))


pdf("CGcontexteDistribution.pdf", width=3.7,height=10)
print(ggarrange(MetTypeod1,MetTypeod2,MetTypeod3,MetTypeovary,ncol=1,nrow=4))
dev.off()

####Turning CpGContext graph --> Stacked Vertical Bar Plot 
###renaming columns so that i can differentiate them in master df 
mergeod1 = flipod1.c %>% rename ("Magnum(od1)"='avg met lvl')
mergeod2 = flipod2.c %>% rename ("Isthmus(od2)"='avg met lvl')
mergeod3 = flipod3.c %>% rename ("Shell Gland (od3)"='avg met lvl')
mergeovary = flipovary.c %>% rename ("Ovary"='avg met lvl')


###making the master df 
mastercgdata <- mergeod1 %>% 
  merge(mergeod2, by = "nucleotide pair type", all = TRUE) %>% 
  merge(mergeod3, by = "nucleotide pair type", all = TRUE) %>%
  merge(mergeovary, by = "nucleotide pair type", all = TRUE)

###isolating CH, CHG, and CHH (isolated rows w/ "nucleotide pair type" values that had the letter H in it )
###mastercgdata2 <- mastercgdata[grepl("H", mastercgdata$`nucleotide pair type`), ]

###isolating CG, CHG, and CHH 
mastercgdata2 <- mastercgdata[mastercgdata$`nucleotide pair type` %in% c("CG", "CHG", "CHH"), ]


###Flipping this df so that it is compatible to be convert from wide-->long farmat (which is req. to plot stacked bar plot)
##FLIP the x & y axis of df 
mastercgdata3 <- t(mastercgdata2)
##for some reason, mastercgdata3 was NOT recognized as a df by R, so I made a new df that recognized the matrix as a df
mastercgdata4 = as.data.frame(mastercgdata3)
##Converting the first row --> column names / header 
mastercgdata5 <- mastercgdata4[-1,]
colnames(mastercgdata5) <- as.character(mastercgdata4[1,])
##converting the rownames to an actual column vector 
mastercgdata6 <- mastercgdata5 %>% 
  rownames_to_column(var = "tissue")

###Converting this df (mastercgdata6) from wide --> long format so that it is easier (compatible) to make stacked bar plot
##converting df from wide-->long format (key = rownames, value = the actual #s in the boxes)
mastercgdata7 = tidyr::gather(mastercgdata6, key= "nucleotide pair type", value = "avg met lvl", -tissue)

# specify the order that x variable will be listed
mastercgdata7$tissue <- factor(mastercgdata7$tissue, levels = c("Magnum(od1)","Isthmus(od2)","Shell Gland (od3)","Ovary"))

##Creating Stacked Barplot 
mean_mCcontexts=ggplot(mastercgdata7, aes(x = tissue, y = as.numeric(`avg met lvl`), fill = `nucleotide pair type`)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55) +
  labs(title = "Methylation Context Contributions Across Tissues",
       x = "Tissue",
       y = "Average Methylation (mean_mC)",
       fill = "Methylation context") + scale_fill_manual(values = c("skyblue","lightgreen","pink")) + theme_bw()+theme(axis.text.x=element_text(angle=30,size=9, vjust=0.88,hjust=0.9),
                                                                                                   axis.text.y=element_text(size=7),
                                                                                                   axis.title.x=element_text(size=11),
                                                                                                   axis.title.y = element_text(size=11),
                                                                                                   plot.margin=margin(0.8,0.8,0.8,2,"cm"),legend.position="right")

#printing the graph 
pdf("MetContext.pdf", width=6,height=6.5)

print(mean_mCcontexts,
      labels=c("Methylation Context Contributions Across Tissues"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()
