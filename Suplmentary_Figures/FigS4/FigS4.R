
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

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/Fig4S")


##################
# FigS4A - MIP coverage plot from vcf files 
##################

##read in the VCF file with vcfR and check the file

FigS4A<- read.vcfR("FigS4A_file", verbose = FALSE)

##vcf Summarization and check metadata

head(FigS4A) 
queryMETA(FigS4A) 
queryMETA(FigS4A, element = 'DP')

head(is.polymorphic(FigS4A, na.omit = TRUE))
head(is.biallelic(FigS4A))
queryMETA(FigS4A, element = 'FORMAT=<ID=DP')
strwrap(FigS4A@meta[1:7])

# If you needed you can Subset samples
FigS4A[,1:10]

# The fix region
# The fix region contains information for each variant which is sometimes summarized over all samples. The first eight columns of the fixed region and are titled CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO. This is per variant information which is ‘fixed’, or the same, over all samples. The first two columns indicate the location of the variant by chromosome and position within that chromosome. Here, the ID field has not been used, so it consists of missing data (NA). The REF and ALT columns indicate the reference and alternate allelic states. When multiple alternate allelic states are present they are delimited with commas. The QUAL column attempts to summarize the quality of each variant over all samples. The FILTER field is not used here but could contain information on whether a variant has passed some form of quality assessment.
tail(getFIX(FigS4A))

# The gt region
#The gt (genotype) region contains information about each variant for each sample. The values for each variant and each sample are colon delimited. Multiple types of data for each genotype may be stored in this manner. The format of the data is specified by the FORMAT column (column nine). Here we see that we have information for GT, AD, DP, GQ and PL. The definition of these acronyms can be referenced by querying the the meta region, as demonstrated previously. Every variant does not necessarily have the same information (e.g., SNPs and indels may be handled differently), so the rows are best treated independently. Different variant callers may include different information in this region.
FigS4A@gt[1:6, 1:4]

#quick check cumulative read depth distribution across samples

DP_samples <- extract.gt(FigS4A, element='DP', as.numeric=TRUE)
rownames(DP_samples ) <- 1:nrow(DP_samples )
head(DP_samples )

heatmap.bp(DP_samples)

is.na(DP_samples [na.omit(DP_samples  == 0)]) <- TRUE

heatmap.bp(log(DP_samples), cbarplot = F, rbarplot = F, min= 0.2, max =0.5, legend = TRUE, clabels = F, rlabels = TRUE, na.rm = TRUE,
           scale = c("column"),
           #col.ramp = viridisLite::viridis(n = 100, alpha = 1))
           #col.ramp = colorRampPalette(c("yellow", "orange", "red"))(100))
           
           col.ramp = colorRampPalette(c("orange", "darkred"))(110))



# Save plot 
ggsave("FigS4A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig4A.pdf", dpi=600, width=7.5, height=7)

##################
# FigS4B - Plot SNP-density from vcf files 
##################

# load the following libraries 
library(qqman)
library(calibrate)
library(CMplot)

# Load SNP density formatted data 

FigS4B<- read.csv("FigS4B_file.csv", header = T) # Look the attached file. This one of intermediate file for selection analysis. It contains four columns Snp_coordinate, Chro_no, SNP_position and P_value (not required for SNP density plot)


#plot SNP-density plot

CMplot(MIP_snp_density,plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)          


