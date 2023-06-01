
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################


####################################################
# Fig4- R script to generate PCA plot using VCF file   
####################################################

##  Please install the R packages first
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

#Load the library
library(gdsfmt)
library(SNPRelate)
library(FactoClass)
library(scatterplot3d)
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(vcfR)
library(tidyverse)
library(dplyr)
library(tidyr)


setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig4")


# select the data

input='Fig4_file.vcf'


#) Please change the name of the prefix of output files if needed
name = "pf_ephi"

# Define Two pops 

pop1_metadata<-read.table("pop1_metadata.csv", sep="\t")
A1<-colnames(pop1_metadata) <- c("sample_id", "hrp23status")


pop12_metadata<-read.table("pop2_metadata.tsv", sep="\t")
A2<-colnames(pop2_metadata) <- c("sample_id", "R622Istatus")


#Please change the colors if needed

cor.hpr23 <- c( "tomato4", "#CC79A7", "#0072B2", "#009E73")

#define shape

shape.R622Istatus <- c(5, 19)

# define vcf file 
vcf.fn <-input

ephi = paste('ephi_', name, '.gds', sep = "")
ephi

# Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, ephi,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(ephi)

#calculate the eigenvectors and eigenvalues for principal component analysis.
ephi_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

View(ephi_pca)
ephi_pca$varprop[1:10]*100


pc.percent <- ephi_pca$varprop*100
head(round(pc.percent, 2))

# eigenvalues
EV=plot(ephi_pca$eigenval[1:10], type = "o", col = "red", xlab = "PCA", ylab = "Eigenvalues",main = "Eigenvalues for EPHI data analysis")

par(mfrow=c(1, 1))
EV
# variance proportion (%)
plot(ephi_pca$varprop[1:10]*100, type = "o", 
        col = "black", siz=1, xlab = "PC", ylab = "(%)",
        main = "EPHI samples Variance Explained")

barplot(ephi_pca$varprop[1:10]*100, type = "o", col = "darkblue", xlab = "PC", ylab = "Percentage of variance explained",
            main = " variance explained ")

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = ephi_pca$eigenvect[,1]    # the first eigenvector
EV2 = ephi_pca$eigenvect[,2]    # the second eigenvector
EV3 = ephi_pca$eigenvect[,3]    # the 3rd eigenvector
EV4 = ephi_pca$eigenvect[,4]    # the 4th eigenvector


## Get a list of pops following the order of the sample IDs and plot pca1 and 2

## Merge pop with you pca file

pop1 = factor(pop1_metadata$hrp23status)[match(ephi_pca$sample.id, 
                                             pop1_metadata$sample_id)]


pop2 = factor(pop2_metadata$R622Istatus)[match(ephi_pca$sample.id, 
                                             pop2_metadata$sample_id)]

# Write pca file.tsv format
tab <- data.frame(ephi_pca$sample.id,pop3, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA1_2hpr23status_final.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")
head(tab)

# plot pca as pdf 

figpdf = paste('Fig4.pdf', sep="")
pdf(file = figpdf)

plot(EV1, EV2,xlab="PC1(7.3%)", ylab="PC2(5.1%)", col=cor.hpr23[as.integer(pop1)],
     pch=shape.R622Istatus[as.integer(pop2)],cex=1)
legend("topleft", legend=levels(pop1), bg="transparent",pch=16, cex=1,col=cor.hpr23,text.col=cor.hpr23)
legend("topright", legend=levels(pop2), bg="transparent",pch=shape.R622Istatus, cex=1)
abline(v=0, h=0, col="black", lwd=1, lty=2)

dev.off()

# plot sgv as pca 

svglite("Fig4.svg", width=7.5, height=7)
plot(EV1, EV2,xlab="PC1(7.3%)", ylab="PC2(5.1%)", col=cor.hpr23[as.integer(pop1)],
     pch=shape.R622Istatus[as.integer(pop2)],cex=1)
legend("topleft", legend=levels(pop1), bg="transparent",pch=16, cex=1,col=cor.hpr23,text.col=cor.hpr23)
legend("topright", legend=levels(pop2), bg="transparent",pch=shape.R622Istatus, cex=1)
abline(v=0, h=0, col="black", lwd=1, lty=2)

dev.off()

sessionInfo()


                                                                                                                        