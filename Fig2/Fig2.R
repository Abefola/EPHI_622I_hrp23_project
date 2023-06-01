
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################


#################
#################
#  Fig2 - CRT_MDR1_622I mutations prevalence 
#################
#################

getwd() # To check the current directory 

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig2")

# Load required libraries 
library (ggplot2)
library(ggpubr)
library(UpSetR)

# Define colors to be used 

bar_cols1 <- c( "black", "#999999","#56B4E9", 
                       "#009E73",  "#E69F00","#0072B2", 
                       "#D55E00", "#CC79A7", "#0072B2")
                       bar_cols2 <- c("#0072B2","darkgreen","grey", "dodgerblue",
                                               "black", "darkcyan", "red","orange","tomato4")
# Load you data   
                       
 Fig2 <- read.csv("Fig2_file.csv", header=TRUE, sep="," )

#Columns should contain
# Sample Id, mutant names, and Wild genotype (if it contains wild alleles for all loci )

#Plot the data 
upset(Fig2, nsets = 10, nintersects = 30, mb.ratio = c(0.5, 0.5),
text.scale = c (2.5, 2.5, 2.5, 2.5, 2, 2),
order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE), 
keep.order = T, main.bar.color = bar_cols1,
point.size = 3.5, line.size = 1)

# Save plot 
 ggsave("Fig2.svg", dpi=600, width=7.5, height=7)
 ggsave("Fig2.pdf", dpi=600, width=7.5, height=7)
                       
                       