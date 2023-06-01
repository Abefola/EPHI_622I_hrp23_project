
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

#################
#################
# Parasites density vs MIP sequencing coverage
#################
#################


getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS2") # change this to your working directory


#Install Library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)

#########
# FigS2A - PCR Pf Parasitaemia Distribution
#########

FigS2A<-read.csv("FigS2A_file.csv", header=T)

### change the Parasitaemia to categorical groups
ggplot(FigS2A, aes(x=pfldh_pcr_parasitemia_mean)) + 
  geom_histogram(aes(y=..density..),log="x",  bindwidth =0.5,outlier.test=F,  colour = "blue", fill = "white") + 
  geom_density( alpha =0.2, fill = "#FF6666")+  labs(x="PCR Parasitaemia", y="Density")+
  ggtitle("PCR Pf Parasitaemia Distribution")

# Save plot 
ggsave("FigS2A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1SA.pdf", dpi=600, width=7.5, height=7)


##############
# FigS2B - Average Barcode Count across Different Parasitaemia Level
##############

FigS2B <- read.csv("FigS2B_file",header=T,  sep=",")

ggplot(FigS2B, aes(x = "", y =log(mean_barcode_count), color=log(pfldh_pcr_parasitemia_mean)))+
  geom_boxplot(outlier.shape=NA)+ # remove outlier from boxplot
  geom_jitter()+
  theme_bw () +
  labs(x="log(Pfldh-PCR Parasitaemia)", y="log( Barcode Count)") +
  scale_color_viridis(option = "D")+
  ggtitle("Average Barcode Count across Different Parasitaemia Level")

# Save plot 
ggsave("FigS2B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1SB2.pdf", dpi=600, width=7.5, height=7)


