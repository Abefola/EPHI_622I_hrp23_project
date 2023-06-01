
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

#################
#################
# FigS7A - crt UpSet plot
#################
#################

getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS7") # change this to your working directory

#Install Library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(UpSetR)

#Load the data 

FigS7A_file <- read.csv("FigS7A_file.csv", header=TRUE, sep="," )

#Columns contain
# Sample Id, mutant names, and Wild genotype (if it contains wild alleles for all loci )

### Define colots 

bar_cols1 <- c( "#56B4E9",
                         "#F0E442", "#CC79A7", "#E69F00",  "darkred")
                         
upset(FigS7A_file, nsets = 9, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE), keep.order = FALSE, main.bar.color = bar_cols1)

# Save plot 
ggsave("FigS7A.svg", dpi=600, width=7.5, height=7)
ggsave("FigS7A.pdf", dpi=600, width=7.5, height=7)



################
################
# FigSBA - crt per hrp23 status
################
################

FigS7B_file<- read.csv('FigS7B_file.csv',header=T,  sep=",")
cbp2 <- c("brown1","chocolate","grey60")

# plot
ggplot(FigS7B,aes(x=Mutation, y=prev, fill = hrp23))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Mutation name", y="Prevalence (%)") +
  scale_fill_discrete(name="hrp23-deletion status") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
scale_fill_manual(values=cbp2)+
ggtitle("crt mutations per hrp23-deletion status")+
theme_bw()

# Save plot 
ggsave("FigS7B.svg", dpi=600, width=7.5, height=7)
ggsave("FigS7B.pdf", dpi=600, width=7.5, height=7)


