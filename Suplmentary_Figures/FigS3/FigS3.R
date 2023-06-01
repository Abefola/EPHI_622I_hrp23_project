
####################################################
# Written by Abebe Fola    
# Date 07/05/2022
####################################################

# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ade4)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(ggplot2)


#) Set your working directory
getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS3")

## read in the VCF file with vcfR and check the file

FigS3<- read.vcfR("FigS3_file.vcf", verbose = FALSE)

## vcf Summarization and check metadata

head(FigS3) 
queryMETA(FigS3) 
queryMETA(FigS3, element = 'DP')

head(is.polymorphic(FigS3, na.omit = TRUE))
head(is.biallelic(FigS3))
queryMETA(FigS3, element = 'FORMAT=<ID=DP')
strwrap(FigS3@meta[1:7])

# If you needed you can Subset samples
FigS3[,1:10]

# The fix region
# The fix region contains information for each variant which is sometimes summarized over all samples. The first eight columns of the fixed region and are titled CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO. This is per variant information which is ‘fixed’, or the same, over all samples. The first two columns indicate the location of the variant by chromosome and position within that chromosome. Here, the ID field has not been used, so it consists of missing data (NA). The REF and ALT columns indicate the reference and alternate allelic states. When multiple alternate allelic states are present they are delimited with commas. The QUAL column attempts to summarize the quality of each variant over all samples. The FILTER field is not used here but could contain information on whether a variant has passed some form of quality assessment.
tail(getFIX(FigS3))

# The gt region
#The gt (genotype) region contains information about each variant for each sample. The values for each variant and each sample are colon delimited. Multiple types of data for each genotype may be stored in this manner. The format of the data is specified by the FORMAT column (column nine). Here we see that we have information for GT, AD, DP, GQ and PL. The definition of these acronyms can be referenced by querying the the meta region, as demonstrated previously. Every variant does not necessarily have the same information (e.g., SNPs and indels may be handled differently), so the rows are best treated independently. Different variant callers may include different information in this region.
FigS3@gt[1:6, 1:4]

#######################
# Calculate missingness 
#############

vcfephi.fn <- "FigS3_file.vcf"

# Reformat (change vcf to gds file )

snpgdsVCF2GDS(vcfephi.fn, "ephisamples.gds", method="biallelic.only")

# Summary
snpgdsSummary("ephisamples.gds")

genofile_ephisamples <- snpgdsOpen("ephisamples.gds")

## get sample id and SNP id

sample.id <- read.gdsn(index.gdsn(genofile_ephisamples, "sample.id"))

snp.id <- read.gdsn(index.gdsn(genofile_ephisamples, "snp.id"))


# calculate sample missingness rate - checks samples containing not calls across the genome

samplemissr <- snpgdsSampMissRate(genofile_ephisamples, sample.id=sample.id, snp.id=snp.id, with.id=FALSE)

#write.csv(samplemissr,"sample_missingness_ephisamples.csv")

figpdf = paste('FigS3A.pdf', sep="")
pdf(file = figpdf)
hist(samplemissr, breaks= 30,xlim= c(0, 1), col="grey",   xlab="Missingness",
     main= "Sample Missingness")

abline(v=0.5, col="red", lwd=3, lty=2)

dev.off()

# calculate SNP missigness rate - checks individual missingness rate across samples
SNPmissrate<-snpgdsSNPRateFreq(genofile_ephisamples, sample.id=sample.id, snp.id=snp.id, with.id=FALSE,
                               with.sample.id=FALSE, with.snp.id=FALSE)

#write.csv(SNPmissrate,"SNP_missingness_ephisamples.csv")

figpdf = paste('FigS3B.pdf', sep="")
pdf(file = figpdf)
hist(SNPmissrate$MissingRate, breaks= 20, col="grey",   xlab="Missingness", main= "SNP Missingness")
abline(v=0.46, col="red", lwd=3, lty=2)

dev.off()


