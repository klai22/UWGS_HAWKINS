##setwd 
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/MetType")

##load packages
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library (biomaRt)

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

##Creating Hyper vs. Hypomethylation plots 
##OD12
#Calc. # of hypermethylated DMRs 
hyperod12=nrow(od1_od2sig[od1_od2sig$t > 0, ])
hypood12=nrow(od1_od2sig[od1_od2sig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod12=c("od1","od2")
methylationtypevalueod12=c( hyperod12, hypood12)
od12hyp = data.frame(methylationtypeod12, methylationtypevalueod12)
#plot the df 
od12plot =ggplot(od12hyp,aes(x=methylationtype,y=methylationtypevalueod12)) + 
  geom_col(aes(fill = methylationtypeod12), width = 0.5) + 
  scale_fill_manual(values = c("pink", "yellow")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od1vs.od2)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")


##OD13
#Calc. # of hypermethylated DMRs 
hyperod13=nrow(od1_od3sig[od1_od3sig$t > 0, ])
hypood13=nrow(od1_od3sig[od1_od3sig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod13=c("od1","od3")
methylationtypevalueod13=c( hyperod13, hypood13)
od13hyp = data.frame(methylationtypeod13, methylationtypevalueod13)
#plot the df 
od13plot = ggplot(od13hyp,aes(x=methylationtypeod13,y=methylationtypevalueod13)) + 
  geom_col(aes(fill = methylationtypeod13), width = 0.5) + 
  scale_fill_manual(values = c("pink", "lightgreen")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od1vs.od3)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill ="Tissue Type")

##OD1o
#Calc. # of hypermethylated DMRs 
hyperod1o=nrow(od1_ovarysig[od1_ovarysig$t > 0, ])
hypood1o=nrow(od1_ovarysig[od1_ovarysig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod1o=c("od1","ovary")
methylationtypevalueod1o=c( hyperod1o, hypood1o)
od1ohyp = data.frame(methylationtypeod1o, methylationtypevalueod1o)
#plot the df 
od1oplot = ggplot(od1ohyp,aes(x=methylationtypeod1o,y=methylationtypevalueod1o)) + 
  geom_col(aes(fill = methylationtypeod1o), width = 0.5) + 
  scale_fill_manual(values = c("pink", "skyblue")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od1vs.ovary)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill ="Tissue Type")


##OD23
#Calc. # of hypermethylated DMRs 
hyperod23=nrow(od2_od3sig[od2_od3sig$t > 0, ])
hypood23=nrow(od2_od3sig[od2_od3sig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod23=c("od2","od3")
methylationtypevalueod23=c( hyperod23, hypood23)
od23hyp = data.frame(methylationtypeod23, methylationtypevalueod23)
#plot the df 
od23plot = ggplot(od23hyp,aes(x=methylationtypeod23,y=methylationtypevalueod23)) + 
  geom_col(aes(fill = methylationtypeod23), width = 0.5) + 
  scale_fill_manual(values = c("yellow", "lightgreen")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Isthmus DMRs (od2vs.od3)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill ="Tissue Type")


##OD2O
#Calc. # of hypermethylated DMRs 
hyperod2o=nrow(od2_ovarysig[od2_ovarysig$t > 0, ])
hypood2o=nrow(od2_ovarysig[od2_ovarysig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod2o=c("od2","ovary")
methylationtypevalueod2o=c( hyperod2o, hypood2o)
od2ohyp = data.frame(methylationtypeod2o, methylationtypevalueod2o)
#plot the df 
od2oplot = ggplot(od2ohyp,aes(x=methylationtypeod2o,y=methylationtypevalueod2o)) + 
  geom_col(aes(fill = methylationtypeod2o), width = 0.5) + 
  scale_fill_manual(values = c("yellow", "skyblue")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Isthmus DMRs (od2vs.ovary)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill ="Tissue Type")
##OD3O
#Calc. # of hypermethylated DMRs 
hyperod3o=nrow(od3_ovarysig[od3_ovarysig$t > 0, ])
hypood3o=nrow(od3_ovarysig[od3_ovarysig$t < 0, ])
#Create df for OD1 graph 
methylationtypeod3o=c("od3","ovary")
methylationtypevalueod3o=c( hyperod3o, hypood3o)
od3ohyp = data.frame(methylationtypeod3o, methylationtypevalueod3o)
#plot the df 
od3oplot = ggplot(od3ohyp,aes(x=methylationtypeod3o,y=methylationtypevalueod3o)) + 
  geom_col(aes(fill = methylationtypeod3o), width = 0.5) + 
  scale_fill_manual(values = c("lightgreen", "skyblue")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Shell Gland DMRs (od3vs.ovary)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill ="Tissue Type")

#combining all graphs???? ERROR???? 
pdf("mettypedata(color coded).pdf", width=20,height=10)
print(ggarrange(od12plot,od13plot,od1oplot,od23plot,od2oplot,od3oplot,
                ncol=3,nrow=2))
dev.off() 

#making inverse duplicate graphs 

##OD21
#Create df for OD21 graph 
methylationtypeod21=c("od2","od1")
methylationtypevalueod21=c(hypood12,hyperod12)
od21hyp = data.frame(methylationtypeod21, methylationtypevalueod21)
# specify the order that x variable will be listed
od21hyp$methylationtypeod21 <- factor(od21hyp$methylationtypeod21, levels = c("od2","od1"))
#plot the df 
od21plot =ggplot(od21hyp,aes(x=methylationtypeod21,y=methylationtypevalueod21)) + 
  geom_col(aes(fill = methylationtypeod21), width = 0.5) + 
  scale_fill_manual(values = c("yellow", "pink")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od2vs.od1)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##OD21
#Create df for OD21 graph 
methylationtypeod21=c("od2","od1")
methylationtypevalueod21=c(hypood12,hyperod12)
od21hyp = data.frame(methylationtypeod21, methylationtypevalueod21)
# specify the order that x variable will be listed
od21hyp$methylationtypeod21 <- factor(od21hyp$methylationtypeod21, levels = c("od2","od1"))
#plot the df 
od21plot =ggplot(od21hyp,aes(x=methylationtypeod21,y=methylationtypevalueod21)) + 
  geom_col(aes(fill = methylationtypeod21), width = 0.5) + 
  scale_fill_manual(values = c("yellow", "pink")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od2vs.od1)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##OD31
#Create df for OD21 graph 
methylationtypeod31=c("od3","od1")
methylationtypevalueod31=c(hypood13,hyperod13)
od31hyp = data.frame(methylationtypeod31, methylationtypevalueod31)
# specify the order that x variable will be listed
od31hyp$methylationtypeod31 <- factor(od31hyp$methylationtypeod31, levels = c("od3","od1"))
#plot the df 
od31plot =ggplot(od31hyp,aes(x=methylationtypeod31,y=methylationtypevalueod31)) + 
  geom_col(aes(fill = methylationtypeod31), width = 0.5) + 
  scale_fill_manual(values = c("lightgreen", "pink")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od3vs.od1)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##OD32
#Create df for OD21 graph 
methylationtypeod32=c("od3","od2")
methylationtypevalueod32=c(hypood23,hyperod23)
od32hyp = data.frame(methylationtypeod32, methylationtypevalueod32)
# specify the order that x variable will be listed
od32hyp$methylationtypeod32 <- factor(od32hyp$methylationtypeod32, levels = c("od3","od2"))
#plot the df 
od32plot =ggplot(od32hyp,aes(x=methylationtypeod32,y=methylationtypevalueod32)) + 
  geom_col(aes(fill = methylationtypeod32), width = 0.5) + 
  scale_fill_manual(values = c("lightgreen", "yellow")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (od3vs.od2)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##ODO1
#Create df for OD21 graph 
methylationtypeodo1=c("ovary","od1")
methylationtypevalueodo1=c(hypood1o,hyperod1o)
odo1hyp = data.frame(methylationtypeodo1, methylationtypevalueodo1)
# specify the order that x variable will be listed
odo1hyp$methylationtypeodo1 <- factor(odo1hyp$methylationtypeodo1, levels = c("ovary","od1"))
#plot the df 
odo1plot =ggplot(odo1hyp,aes(x=methylationtypeodo1,y=methylationtypevalueodo1)) + 
  geom_col(aes(fill = methylationtypeodo1), width = 0.5) + 
  scale_fill_manual(values = c("skyblue", "pink")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (ovary vs.od1)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##ODO2
#Create df for OD21 graph 
methylationtypeodo2=c("ovary","od2")
methylationtypevalueodo2=c(hypood2o,hyperod2o)
odo2hyp = data.frame(methylationtypeodo2,methylationtypevalueodo2)
# specify the order that x variable will be listed
odo2hyp$methylationtypeodo2 <- factor(odo2hyp$methylationtypeodo2, levels = c("ovary","od2"))
#plot the df 
odo2plot =ggplot(odo2hyp,aes(x=methylationtypeodo2,y=methylationtypevalueodo2)) + 
  geom_col(aes(fill = methylationtypeodo2), width = 0.5) + 
  scale_fill_manual(values = c("skyblue", "yellow")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (ovary vs.od2)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")

##ODO3
#Create df for OD21 graph 
methylationtypeodo3=c("ovary","od3")
methylationtypevalueodo3=c(hypood3o,hyperod3o)
odo3hyp = data.frame(methylationtypeodo3,methylationtypevalueodo3)
# specify the order that x variable will be listed
odo3hyp$methylationtypeodo3 <- factor(odo3hyp$methylationtypeodo3, levels = c("ovary","od3"))
#plot the df 
odo3plot =ggplot(odo3hyp,aes(x=methylationtypeodo3,y=methylationtypevalueodo3)) + 
  geom_col(aes(fill = methylationtypeodo3), width = 0.5) + 
  scale_fill_manual(values = c("skyblue", "lightgreen")) + 
  xlab("Hypermethylated in...")+
  ylab("DMRs detected") +
  theme_classic() +
  ggtitle("Methylation Levels in Magnum DMRs (ovary vs.od3)")+
  theme(plot.title=element_text(face="italic",hjust=0.5,size=10),
        axis.text.x=element_text(angle=55,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  labs(fill = "Tissue Type")
#printing out the sets of graphs by tissue 
pdf("mettypedata(inverseduplicates).pdf", width=20,height=22)
print(ggarrange(od12plot,od13plot,od1oplot,od21plot,od23plot,od2oplot,od31plot,od32plot,od3oplot,odo1plot,odo2plot,odo3plot,
                ncol=3,nrow=4))
dev.off() 



#showing all values 
od12hyp
od13hyp
od1ohyp
od23hyp
od2ohyp
od3ohyp




#Making New Graph to show total DMR Count  / Sample 
od12totalcount=nrow(od1_od2sig)
od13totalcount=nrow(od1_od3sig)
od1ototalcount=nrow(od1_ovarysig)
od23totalcount=nrow(od2_od3sig)
od2ototalcount=nrow(od2_ovarysig)
od3ototalcount=nrow(od3_ovarysig)
  
#Creating a new DF w/ each of these totalcount values 
totaldmrcounts <- data.frame(
  tissuescompared = c("od12","od13", "od1o", "od23", "od2o", "od3o"),
  totalDMRcount = c(od12totalcount, od13totalcount, od1ototalcount, od23totalcount, od2ototalcount, od3ototalcount
  ))

#Creating bar graph 
totalDMRcountplot=ggplot(totaldmrcounts,aes(x=tissuescompared,y=totalDMRcount)) + 
  geom_col(fill="lightgreen",width=0.5) + 
  xlab("Tissues Compared") + 
  ylab("Total DMRs detected") +
  theme_classic() +
  ggtitle("Total DMRs detected per Pairwise Comparison") +
  theme(plot.title=element_text(face="italic",hjust=0.5,size=12),
        axis.text.x=element_text(angle=18,size=9, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.margin=margin(1,2,1,2,"cm"))+ 
  scale_x_discrete(labels=c("od12"="od1 vs. od2",
                            "od13"="od1 vs. od3",
                            "od1o"="od1 vs. ovary",
                            "od23"="od2 vs. od3",
                            "od2o"="od2 vs. ovary",
                            "od3o"="od3 vs. od3"))+
  geom_text(aes(label = totalDMRcount), vjust = -0.5, size = 3)

#printing bar graph 
pdf("totalDMRcounts.pdf")
print(totalDMRcountplot,ncol=1,nrow=1)
dev.off() 

