#setting working directory
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/CpGIslands")

#loading packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")

#import data 
od12islands = read.csv("od12CpGislands.csv", header=FALSE)
od12shores = read.csv("od12CpGshores.csv", header=FALSE)
od12shelves = read.csv("od12CpGshelves.csv", header=FALSE)
od13islands = read.csv("od13CpGislands.csv", header=FALSE)
od13shores = read.csv("od13CpGshores.csv", header=FALSE)
od13shelves = read.csv("od13CpGshelves.csv", header=FALSE)
od1oislands = read.csv("od1oCpGislands.csv", header=FALSE)
od1oshores = read.csv("od1oCpGshores.csv", header=FALSE)
od1oshelves = read.csv("od1oCpGshelves.csv", header=FALSE)
od23islands = read.csv("od23CpGislands.csv", header=FALSE)
od23shores = read.csv("od23CpGshores.csv", header=FALSE)
od23shelves = read.csv("od23CpGshelves.csv", header=FALSE)
od2oislands = read.csv("od2oCpGislands.csv", header=FALSE)
od2oshores = read.csv("od2oCpGshores.csv", header=FALSE)
od2oshelves = read.csv("od2oCpGshelves.csv", header=FALSE)
od3oislands = read.csv("od3oCpGislands.csv", header=FALSE)
od3oshores = read.csv("od3oCpGshores.csv", header=FALSE)
od3oshelves = read.csv("od3oCpGshelves.csv", header=FALSE)

od12totalDMRs = read.table("od12CGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od13totalDMRs = read.table("od13CGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od1ototalDMRs = read.table("od1oCGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od23totalDMRs = read.table("od23CGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od2ototalDMRs = read.table("od2oCGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od3ototalDMRs = read.table("od3oCGisland.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#Create DF 
#row counting (creating %s  out of the counts )
od12islandscount = nrow(od12islands) / nrow(od12totalDMRs) * 100
od12shorescount = nrow(od12shores) / nrow(od12totalDMRs) * 100
od12shelvescount = nrow(od12shelves) / nrow(od12totalDMRs) * 100
od13islandscount = nrow(od13islands) / nrow(od13totalDMRs) * 100
od13shorescount = nrow(od13shores) / nrow(od13totalDMRs) * 100
od13shelvescount = nrow(od13shelves) / nrow(od13totalDMRs) * 100
od1oislandscount = nrow(od1oislands) / nrow(od1ototalDMRs) * 100
od1oshorescount = nrow(od1oshores) / nrow(od1ototalDMRs) * 100
od1oshelvescount = nrow(od1oshelves) / nrow(od1ototalDMRs) * 100
od23islandscount = nrow(od23islands) / nrow(od23totalDMRs) * 100
od23shorescount = nrow(od23shores) / nrow(od23totalDMRs) * 100
od23shelvescount = nrow(od23shelves) / nrow(od23totalDMRs) * 100
od2oislandscount = nrow(od2oislands) / nrow(od2ototalDMRs) * 100
od2oshorescount = nrow(od2oshores) / nrow(od2ototalDMRs) * 100
od2oshelvescount = nrow(od2oshelves) / nrow(od2ototalDMRs) * 100
od3oislandscount = nrow(od3oislands) / nrow(od3ototalDMRs) * 100
od3oshorescount = nrow(od3oshores) / nrow(od3ototalDMRs) * 100
od3oshelvescount = nrow(od3oshelves) / nrow(od3ototalDMRs) * 100

#Creating the DF 
cpgislanddata <- data.frame(
  tissuescompared = c("od12", "od12","od12","od13", "od13", "od13", "od1o", "od1o", "od1o", "od23", "od23", "od23", "od2o", "od2o", "od2o", "od3o","od3o", "od3o"),
  cpgcontext = c("islands","shores","shelves","islands","shores","shelves","islands","shores","shelves","islands","shores","shelves","islands","shores","shelves","islands","shores","shelves"),
  DMRcount = c(od12islandscount, od12shorescount, od12shelvescount, od13islandscount, od13shorescount, od13shelvescount, od1oislandscount, od1oshorescount, od1oshelvescount, od23islandscount, od23shorescount, od23shelvescount, od2oislandscount, od2oshorescount, od2oshelvescount, od3oislandscount, od3oshorescount, od3oshelvescount
))


#Plotting Bargraph 
cpgislandplot=ggplot(cpgislanddata, aes(x = tissuescompared, y = DMRcount, fill = cpgcontext)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightgreen","skyblue","pink"), labels=c("islands"="cpgislands",
                                                                        "shelves"="cpgshelves",
                                                                        "shores"="cpgshores")) +
  labs(x = "Tissues Compared", y = "DMR Count (% of total DMRs)", fill = "CpGIsland Context")+theme_classic()+ggtitle("DMR CpG Island Context Distribution")+
  theme(axis.text.x=element_text(angle=25,size=8, vjust=0.88,hjust=0.9),
        axis.text.y=element_text(size=7),
        axis.title.x=element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.margin=margin(0.8,0.8,0.8,2,"cm"),legend.position="right")+ 
  scale_x_discrete(labels=c("od12"="od1 vs. od2",
                             "od13"="od1 vs. od3",
                             "od1o"="od1 vs. ovary",
                             "od23"="od2 vs. od3",
                             "od2o"="od2 vs. ovary",
                             "od3o"="od3 vs. ovary"))
#printing plot 

pdf("CpGislands(percent).pdf", width=10,height=8)

print(cpgislandplot,
      labels=c("DMR CpG Island Context Distribution"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()

