####################################################
# Written by Abebe Fola    
# Date 07/05/2022
####################################################

##################
# FigS5A - COI analysis
##################

# Load library 
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(dplyr)
library(tidyr)
library(viridis)
library(leaflet)
library(rhandsontable)
library(sp)
library(rgeos)
library(ggpubr)
library(scales)
library(gridExtra)
library(scatterplot3d)
library(knitr)
library(ggExtra)
library(GGally)
library(viridis)
library("ggsci")
library(moimix)
library(flexmix)
library(lattice)
library("ggpubr")

getwd()
setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS5")

FigS5A<-read.csv("FigS5A_file", header=T,  sep=",")

# COI
ggplot(FigS5A, aes(coi_mean)) + geom_histogram(binwidth = 0.5, fill = "gray43")+
  labs(y="Number of Samples", x="COI") 


# Save plot 
ggsave("FigS5A.svg", dpi=600, width=7.5, height=7)
ggsave("FigS5A.pdf", dpi=600, width=7.5, height=7)

#################
# FigS5B FWS analysis
###############

FigS5B='FigS5B_file.vcf'

vcf_header <- seqVCF_Header(FigS5B)
vcf_header
# recode header format for AD 
vcf_header$format$Number[vcf_header$ID == "AD"] <- "."

# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "ANN")

# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

GDS_EPHI_samples<-seqVCF2GDS(FigS5B,
                             "EPH_samples.gds",
                             header=vcf_header, info.import=info.import,
                             fmt.import=format.import)

seqSummary(GDS_EPHI_samples)

EPHI_isolates <- seqOpen(GDS_EPHI_samples)
seqSummary(EPHI_isolates)


#save sample identifiers

sample.id <- seqGetData(EPHI_isolates, "sample.id")

summary(sample.id)

list (sample.id)

coords <- getCoordinates(EPHI_isolates)

head(coords)

#
set.seed(3334568) 
#An alternative way of assessing MOI is to estimate the FWS statistic. Fws<0.95 is indicative of multi-clonal infection- Refe(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738909/)

fws<- getFws(EPHI_isolates)

#write.csv(FigS5C, "fws_file.csv") # add province info on this file and using pop file above and use mod file for plotting

fws<-read.csv('fws_file.csv', header=T,  sep=",") # remove three samples (Central, Lusaka and Southern province) and order them per province. 

#plot fws
figpdf = paste('FigS5B.pdf', sep="")
pdf(file = figpdf)
boxplot(EPHI_fws$fws, col = "grey", ylab= "FWS")
dev.off()


#################
# FigS5C Mean COI per district 
###############

FigS5C<-read.csv("FigS5C_file.csv", header=T,  sep=",")

dbs2 <- FigS5C %>%
  filter(!is.na(District)) %>%
  group_by(as.factor(District)) %>%
  summarize(av = mean(MOI),
            sd = sd(value),
            n = n()) %>%
  #elements = paste0("n=",n_distinct(comparison))) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = av - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = av + qt(1 - (0.05 / 2), n - 1) * se)

ggplot(dbs2,aes(as.factor(`as.factor(District)`),av)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.3,
                position=position_dodge(0.05))  +
  ylab("Mean COI") + xlab("District") +
  # geom_text(aes(y=0.03,label=elements),vjust=0.5,size=3) +
  #scale_x_discrete(labels=c("0","1","2","3","4","5")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot 
ggsave("FigS5C.svg", dpi=600, width=7.5, height=7)
ggsave("FigS5C.pdf", dpi=600, width=7.5, height=7)


