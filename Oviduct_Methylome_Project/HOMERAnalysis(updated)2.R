#setting working directory
setwd("/Users/kennethlai/desktop/HAWKINS/USDA_Project/Analysis/UPDATED_PARAMETERS/Homer")

#loading packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")

#import data 
od12 = read.csv("od12DMRcontexts.csv")
od13 = read.csv("od13DMRcontexts.csv")
od1o = read.csv("od1oDMRcontexts.csv")
od23 = read.csv("od23DMRcontexts.csv")
od2o = read.csv("od2oDMRcontexts.csv")
od3o = read.csv("od3oDMRcontexts.csv")

od12totalDMRs = read.table("od12homerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od13totalDMRs = read.table("od13homerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od1ototalDMRs = read.table("od1ohomerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od23totalDMRs = read.table("od23homerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od2ototalDMRs = read.table("od2ohomerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
od3ototalDMRs = read.table("od3ohomerupdated.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#looking at all the possible values found in the Annotation Category 
unique (od12$Annotation)

#counting the # of rows associated w/ each category
##OD12
od12introncount <- sum(grepl("^intron", od12$Annotation)) / nrow(od12totalDMRs) * 100
cat("# of introns: ", od12introncount)

od12exoncount <- sum(grepl("^exon", od12$Annotation)) / nrow(od12totalDMRs) * 100
cat("# of exons: ", od12exoncount)

od12intergeniccount <- sum(grepl("^Intergenic", od12$Annotation)) / nrow(od12totalDMRs) * 100
cat("# of intergenic: ", od12intergeniccount)

od12promotercount <- sum(grepl("^promoter", od12$Annotation)) / nrow(od12totalDMRs) * 100
cat("# of promoters: ", od12promotercount)

od12TTScount <- sum(grepl("^TTS", od12$Annotation)) / nrow(od12totalDMRs) * 100
cat("# of TTS: ", od12TTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod12= od12introncount + od12exoncount + od12intergeniccount+od12promotercount+od12TTScount
print(sumod12)
nrow(od12)

##OD13
od13introncount <- sum(grepl("^intron", od13$Annotation)) / nrow(od13totalDMRs) * 100
cat("# of introns: ", od13introncount)

od13exoncount <- sum(grepl("^exon", od13$Annotation)) / nrow(od13totalDMRs) * 100
cat("# of exons: ", od13exoncount)

od13intergeniccount <- sum(grepl("^Intergenic", od13$Annotation)) / nrow(od13totalDMRs) * 100
cat("# of intergenic: ", od13intergeniccount)

od13promotercount <- sum(grepl("^promoter", od13$Annotation)) / nrow(od13totalDMRs) * 100
cat("# of promoters: ", od13promotercount)

od13TTScount <- sum(grepl("^TTS", od13$Annotation)) / nrow(od13totalDMRs) * 100
cat("# of TTS: ", od13TTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod13= od13introncount + od13exoncount + od13intergeniccount+od13promotercount+od13TTScount
print(sumod13)
nrow(od13)


#OD1O 
od1ointroncount <- sum(grepl("^intron", od1o$Annotation)) / nrow(od1ototalDMRs) * 100
cat("# of introns: ", od1ointroncount)

od1oexoncount <- sum(grepl("^exon", od1o$Annotation)) / nrow(od1ototalDMRs) * 100
cat("# of exons: ", od1oexoncount)

od1ointergeniccount <- sum(grepl("^Intergenic", od1o$Annotation)) / nrow(od1ototalDMRs) * 100
cat("# of intergenic: ", od1ointergeniccount)

od1opromotercount <- sum(grepl("^promoter", od1o$Annotation)) / nrow(od1ototalDMRs) * 100
cat("# of promoters: ", od1opromotercount)

od1oTTScount <- sum(grepl("^TTS", od1o$Annotation)) / nrow(od1ototalDMRs) * 100
cat("# of TTS: ", od1oTTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod1o= od1ointroncount + od1oexoncount + od1ointergeniccount+od1opromotercount+od1oTTScount
print(sumod1o)
nrow(od1o)


#OD23
od23introncount <- sum(grepl("^intron", od23$Annotation)) / nrow(od23totalDMRs) * 100
cat("# of introns: ", od23introncount)

od23exoncount <- sum(grepl("^exon", od23$Annotation)) / nrow(od23totalDMRs) * 100
cat("# of exons: ", od23exoncount)

od23intergeniccount <- sum(grepl("^Intergenic", od23$Annotation)) / nrow(od23totalDMRs) * 100
cat("# of intergenic: ", od23intergeniccount)

od23promotercount <- sum(grepl("^promoter", od23$Annotation)) / nrow(od23totalDMRs) * 100
cat("# of promoters: ", od23promotercount)

od23TTScount <- sum(grepl("^TTS", od23$Annotation)) / nrow(od23totalDMRs) * 100
cat("# of TTS: ", od23TTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod23= od23introncount + od23exoncount + od23intergeniccount+od23promotercount+od23TTScount
print(sumod23)
nrow(od23)


#OD2O 
od2ointroncount <- sum(grepl("^intron", od2o$Annotation)) / nrow(od2ototalDMRs) * 100
cat("# of introns: ", od2ointroncount)

od2oexoncount <- sum(grepl("^exon", od2o$Annotation)) / nrow(od2ototalDMRs) * 100
cat("# of exons: ", od2oexoncount)

od2ointergeniccount <- sum(grepl("^Intergenic", od2o$Annotation)) / nrow(od2ototalDMRs) * 100
cat("# of intergenic: ", od2ointergeniccount)

od2opromotercount <- sum(grepl("^promoter", od2o$Annotation)) / nrow(od2ototalDMRs) * 100
cat("# of promoters: ", od2opromotercount)

od2oTTScount <- sum(grepl("^TTS", od2o$Annotation)) / nrow(od2ototalDMRs) * 100
cat("# of TTS: ", od2oTTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod2o= od2ointroncount + od2oexoncount + od2ointergeniccount+od2opromotercount+od2oTTScount
print(sumod2o)
nrow(od2o)


#OD3O 
od3ointroncount <- sum(grepl("^intron", od3o$Annotation)) / nrow(od3ototalDMRs) * 100
cat("# of introns: ", od3ointroncount)

od3oexoncount <- sum(grepl("^exon", od3o$Annotation)) / nrow(od3ototalDMRs) * 100
cat("# of exons: ", od3oexoncount)

od3ointergeniccount <- sum(grepl("^Intergenic", od3o$Annotation)) / nrow(od3ototalDMRs) * 100
cat("# of intergenic: ", od3ointergeniccount)

od3opromotercount <- sum(grepl("^promoter", od3o$Annotation)) / nrow(od3ototalDMRs) * 100
cat("# of promoters: ", od3opromotercount)

od3oTTScount <- sum(grepl("^TTS", od3o$Annotation)) / nrow(od3ototalDMRs) * 100
cat("# of TTS: ", od3oTTScount)

#Calc. Sum of categories = 10057, make sure it matches up w/ nrow = 10057
sumod3o= od3ointroncount + od3oexoncount + od3ointergeniccount+od3opromotercount+od3oTTScount
print(sumod3o)
nrow(od3o)

#Creating the DF 
contextdata <- data.frame(
  tissuescompared = c("od12", "od12","od12","od12","od12","od13", "od13", "od13", "od13", "od13", "od1o", "od1o", "od1o", "od1o", "od1o", "od23", "od23", "od23", "od23", "od23", "od2o", "od2o", "od2o", "od2o", "od2o", "od3o","od3o", "od3o", "od3o", "od3o"),
  context = c("Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site", "Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site", "Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site", "Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site", "Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site", "Exon", "Intron", "Intergenic Region", "Promoter/TSS","Transcription Termination Site"),
  DMRcount = c(od12exoncount, od12introncount, od12intergeniccount, od12promotercount, od12TTScount, od13exoncount, od13introncount, od13intergeniccount, od13promotercount, od13TTScount, od1oexoncount, od1ointroncount, od1ointergeniccount, od1opromotercount, od1oTTScount, od23exoncount, od23introncount, od23intergeniccount, od23promotercount, od23TTScount, od2oexoncount, od2ointroncount, od2ointergeniccount, od2opromotercount, od2oTTScount, od3oexoncount, od3ointroncount, od3ointergeniccount, od3opromotercount, od3oTTScount )
)


#Plotting Bargraph 
DMRplot=ggplot(contextdata, aes(x = tissuescompared, y = DMRcount, fill = context)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("purple","orange","lightgreen","skyblue","pink")) +
  labs(x = "Tissues Compared", y = "DMR Count (% of total DMRs)", fill = "DMR Context")+theme_classic()+ggtitle("DMR Genomic Context Distribution")+
  theme(axis.text.x=element_text(angle=25,size=9, vjust=0.88,hjust=0.9),
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

pdf("DMRContext(percent).pdf", width=10,height=8)

print(DMRplot,
      labels=c("DMR Genomic Context Distribution"),
      ncol=1,nrow=1,legend = NULL,
      common.legend = TRUE,
      legend.grob = NULL)

dev.off()











