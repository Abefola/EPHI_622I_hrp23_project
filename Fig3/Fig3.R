

#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################


#################
#################
#  Fig3A - comparison  622I prev per hrp23-deleted vs null parasites
#################
#################

getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig3")

# Load required libraries 
library (ggplot2)
library(ggpubr)
library(UpSetR)
library(lmPerm)
library(coin)
library(gtools)

#Load data 
Fig3A <- read.csv("Fig3A_file.csv")

Fig3A$group <- factor(Fig3A$group)
print(R622I)

# compare mean prevalence difference hrp23-deletion-status
my_comparisons <- list( c("pfhrp2-/3-","pfhrp2+/3+"))


print(ggplot(Fig3A,aes(group,Within_group_k13_622I_prevalence_per_District
))
+ geom_boxplot(fill="grey")
+ stat_sum(alpha=0.7)
# theme_bw()
+ ggtitle("R622I prevalence per Pfhrp23 status at District level")
+ labs(x="Pfhrp23 Status", y="R622I Prevalence (%)") 
+ scale_size(breaks=1:3, range=c(4,8))
+ stat_compare_means(method ="t.test", exact = FALSE )
)

# Save plot 
ggsave("Fig3A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig3A.pdf", dpi=600, width=7.5, height=7)


#################
#################
# Fig3B correlation between 622I vs hrp23-deletions at district level 
# Get your data - Weighted prevalence of 622I and hrp23-deletions 
#################
#################

Fig3B <- read.csv('Fig3B_file.csv', header=T,  sep=",")

attach(Fig3B)


cor2 <- ggscatter(Fig3B, x = "R622I_prev", y = "hrp23_prev", size=5,col="co622Ihrp23",
                  #add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", size=20, fill = "gray"), # Customize reg. line
                  conf.int = TRUE,  label = "District", repel = F,
                  check_overlap = TRUE, hjust = 0, nudge_y=5, 
                  nudge_x = -0.5, xlim=c(0,35),
                  xlab = "Pfhrp2-/3- Frequency (%)", 
                  ylab = "K13 622I Frequency (%))" # Add confidence interval
) +
  
ggtitle("Correlation between 622I and pfhrp2-/3-deletion")
# Add correlation coefficient
cor2 + stat_cor(method = "spearman", label.x = 1, label.y = 32, cex=5) + 
  scale_color_manual(values=c("#1B9E77", "#D95F02"), name= "Co-occurance") 


# Save plot 
ggsave("Fig3B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig3B.pdf", dpi=600, width=7.5, height=7)


# Combine figures
library(pdftools)

par(mfrow = c(1, 2))

plt<-pdftools::pdf_combine(input =
                        list.files(full.names=TRUE,pattern=".pdf"),
                      output = "Fig3.pdf")
