#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

####################################################
# R script to generate PCA plot using VCF file
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

# Set working directory 

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS8")


# load the data

input='FigS8and9_file.vcf'


#) Please change the name of the prefix of output files if needed
name = "pf_ephi"

#) Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name

# Diffine Two pops

ephipopreg<-read.table("ephipopreg_file.tsv", sep="\t")
A1<-colnames(ephipopreg) <- c("sample_id", "Region")


ephipop622I<-read.table("ephipop622I.tsv", sep="\t")
A2<-colnames(ephipop622I) <- c("sample_id", "R622Istatus")


#Please change the colors if needed

cor.reg <- c("#1B9E77", "#D95F02", "#7570B3")

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

#View(ephi_pca)
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

########
########
# FigS8A
########
########

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = ephi_pca$eigenvect[,1]    # the first eigenvector
EV2 = ephi_pca$eigenvect[,2]    # the second eigenvector
EV3 = ephi_pca$eigenvect[,3]    # the 3rd eigenvector
EV4 = ephi_pca$eigenvect[,4]    # the 4th eigenvector


## Get a list of pops following the order of the sample IDs and plot pca 1 and 2

## Merge you pop with you pca file

pop1 = factor(ephipopreg$Region)[match(ephi_pca$sample.id, 
                                       ephipopreg$sample_id)]

pop2 = factor(ephipop622I$R622Istatus)[match(ephi_pca$sample.id, 
                                             ephipop622I$sample_id)]
#head (pop2)
# Write pca file.tsv format
tab <- data.frame(ephi_pca$sample.id,pop1, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA1_2Region.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")
head(tab)

# plot pca as pdf 

figpdf = paste('Fig8SA.pdf', sep="")
pdf(file = figpdf)

plot(EV1, EV2,xlab="PC1(7.3%)", ylab="PC2(5.1%)", col=cor.reg[as.integer(pop1)],
     pch=shape.R622Istatus[as.integer(pop2)],cex=1)
legend("topleft", legend=levels(pop1), bg="transparent",pch=16, cex=1,col=cor.reg,text.col=cor.reg)
legend("topright", legend=levels(pop2), bg="transparent",pch=shape.R622Istatus, cex=1)
abline(v=0, h=0, col="black", lwd=1, lty=2)

dev.off()

########
########
# FigS8B
########
########

barplot(ephi_pca$varprop[1:10]*100, type = "o", col = "darkblue", xlab = "PC", ylab = "Percentage of variance explained",
        main = " variance explained ")

# Save plot 
ggsave("FigS8B.svg", dpi=600, width=7.5, height=7)
ggsave("FigS8B.pdf", dpi=600, width=7.5, height=7)

########
########
# FigS9
########
########

SnpLoad <- snpgdsPCASNPLoading(ephi_pca, genofile)

########
# FigS9A - PC1 SnpLoad value
########

plot(SnpLoad$snploading[1,], col="lightblue", cex= 0.7, pch=19, xlab= "SNPs", 
     ylab="PC1 Loadings",main = "Contribution of each SNP EPHI samples clustering"  )
abline(h= -0.04, col="red", cex= 0.5, lty = "dashed")
abline(h= 0.0, col="black", cex= 0.5, lty = "dashed")
abline(h= 0.04, col="red", cex= 0.5, lty = "dashed")


########
# FigS9B - PC2 SnpLoad value
########

plot(SnpLoad$snploading[2,], col="lightblue", cex= 0.7, pch=19, xlab= "SNPs", 
     ylab="PC2 Loadings",main = "Contribution of each SNP EPHI samples clustering"  )
abline(h= -0.05, col="red", cex= 0.5, lty = "dashed")
abline(h= 0.0, col="black", cex= 0.5, lty = "dashed")
abline(h= 0.05, col="red", cex= 0.5, lty = "dashed")


