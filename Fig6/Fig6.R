

#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

####################################################
# Fig6- IBD analysis   
####################################################

getwd()
setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig6")


deps <- c("tidyverse", "vcfR", "MIPanalyzer", "hmmibdr", "sf", "cowplot", "tidygraph", "ggraph")

# Load library

library(raster)
library(vcfR)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(MIPanalyzer)
library(lubridate)
library(splines2)
library(babynames)
library(sp)
library(sf)
library(ggmap)
library(ggspatial)
library(maps)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(maptools)
library(rgdal)
library(ggthemes)
library(igraph)
library(viridisLite)
library(streamgraph)
library(rhandsontable)
library(stats)
library(rlang)
library(magrittr)
library(ggpubr)
source("R/utils.R")

#######################
#######################
# estimate IBD 
#######################
#######################
# Read in vcf file with only choma samples 
invisible(vcf <- read.vcfR("Only_COI1_EPHI_target_IBC_biallelicsnpmiss05sample05.recode.vcf",verbose=FALSE))

test <- vcf2mipanalyzer_biallelic(vcfR=vcf,verbose=FALSE)
invisible(inf <- inbreeding_mle(test, f=seq(0,1,l=100), report_progress=FALSE,ignore_het=TRUE))

# transofrm ibd into a data frame (var1 and var2 are the samples being compared, but referring to the index of their order in the vcf file)
inf.m <- inf$mle
inf.m2 <- melt(inf.m)
inf.m3 <- inf.m2 %>%
  drop_na()

#get order of samples in vcf file
sampl1 <- test$samples
sampl1 <- sampl1 %>%
  add_rownames()
colnames(sampl1) <- c("Var1","seqID")

sampl2 <- sampl1
colnames(sampl2) <- c("Var2","seqID2")

#merge Var 1 ibd estimates with sample names, by index 
inf.m4 <- merge(sampl1,inf.m3,by=c("Var1"))
#merge Var 2 ibd estimates with sample names, by index 
inf.m5 <- merge(inf.m4,sampl2,by=c("Var2"))


#######################
#Fig6A. IBD comparison per hrp23-deletion status
#######################

Fig6A <- read.csv("Fig6A_file.csv") 


colors=c("#009E73", "#0072B2", "#CC79A7", "tomato4")

my_comparisons1 <- list( c("pfhrp2+/3+","pfhrp2+/3-"), c("pfhrp2+/3+","pfhrp2-/3+"),c("pfhrp2+/3+","pfhrp2-/3-"),  
                         c ("pfhrp2-/3+","pfhrp2+/3-"), c("pfhrp2-/3+","pfhrp2-/3-"),c("pfhrp2+/3-","pfhrp2-/3-"))

ggbox1<- ggboxplot(Fig6A, x = "hrp23_status_f", y = "value", fill= "hrp23_status_f", cex=0.5)
ggbox1+
  scale_y_continuous()+
  theme_gray()+
  scale_fill_manual(values=c("#009E73", "#0072B2", "#CC79A7", "tomato4"))+
  #theme_bw () +
  labs(x="Pfhrp23 status", y="Pairwise IBD value")  +
  
  #theme(axis.text.x = element_text(angle =90))+ 
  stat_compare_means(comparisons=my_comparisons1 )+ 
  # Default method = "kruskal.test" for multiple groups # Global p-value
  ggtitle("IBD value per Pfhrp23 deletion across Ethiopia")

# Save plot 
ggsave("Fig6A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig6A.pdf", dpi=600, width=7.5, height=7)


#######################
#Fig6B. Networks for clonal parasites per hrp23 status
#######################

Fig6B<-read.csv("Fig6B_file.csv", header=TRUE)

ibd4 <- Fig6B

ibd5 <- ibd4 %>%
  filter(value ==1)
relations <- ibd5 %>%
  
  dplyr::dplyr::select(seqID,seqID2,value,hrp23_status_f)

for.checks <- ibd5 %>%
  dplyr::dplyr::select(seqID,seqID2,value,hrp23_status_f,hrp23_status_f2)

relations <- ibd5 %>%
  dplyr::dplyr::select(seqID,seqID2,value,hrp23_status_f)

colnames(relations) <- c("from","to","F","hrp23_status_f")

relations <- relations[order(relations$hrp23_status_f),]
relations$hrp23_status_f <- as.factor(relations$hrp23_status_f)


ibd5 <- ibd4 %>%
  filter(value==1 )
relations <- ibd5 %>%
  dplyr::dplyr::select(seqID,seqID2,value,hrp23_status_f)

colnames(relations) <- c("target","source","F","hrp23_status_f")
relations <- relations[order(relations$hrp23_status_f),]

for.checks <- ibd5 %>%
  dplyr::dplyr::select(seqID,seqID2,value,hrp23_status_f,hrp23_status_f2)
colnames(for.checks) <- c("target","source","F","hrp23_status_f","hrp23_status_f2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links <- relations %>% dplyr::dplyr::select(target, source, F)
#write.csv(links, "linksibd1coi1.csv" )

seq1 <- for.checks %>% dplyr::dplyr::select(target,hrp23_status_f) 
seq2 <- for.checks %>% dplyr::dplyr::select(source,hrp23_status_f2) %>% rename(target=source,hrp23_status_f=hrp23_status_f2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,hrp23_status_f)
nodes <- nodes[order(nodes$hrp23_status_f),]
nodes <- droplevels(nodes)

#ta100mod <- metadcomplete %>% 
 # filter(metadcomplete$seqID %in% nodes$target) %>%
 # dplyr::dplyr::select(seqID, hrp23_status_f)


#write.csv(ta100mod , "ta_ibd1mod.csv") # modify and order per pop
ta100m<-read.csv("ta_ibd1mod.csv", header = T)


gm <- graph_from_data_frame(d=links, vertices=ta100m, directed=F) 


colors <- c(rep("tomato4",68), rep("#CC79A7",5),rep("#0072B2",67),rep("#009E73", 48))

plot(gm, vertex.label=NA,vertex.size=9, main="IBD>=0.95, n=185", vertex.color=colors)
legend(x="topleft", legend=c("hrp2-/3-","hrp2-/3+", "hrp2+/3-", "hrp2+/3+"),
       col=c( "tomato4", "#CC79A7", "#0072B2", "#009E73"), cex=1.2, pch=c(19))

# Save plot 
ggsave("Fig6B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig6B.pdf", dpi=600, width=7.5, height=7)

#######################
#Fig6C - Networks for pfhrp2-/3- clonal parasites per District hrp23_status_f == "pfhrp2-/3-"
#######################

getwd()
Fig6C <- read.csv("Fig6C_file.csv")
ibd2 <- Fig6C
 
ibd3 <- ibd2 %>%
  filter(value == 1.0) %>%
filter (hrp23_status_f == "pfhrp2-/3-")
relations <- ibd3 %>%
 
  dplyr::select(seqID,seqID2,value,District)

for.checks <- ibd3 %>%
  dplyr::select(seqID,seqID2,value,District,District2)


  ibd3 <- ibd2 %>%
  filter(value == 1.0) %>%
filter (hrp23_status_f == "pfhrp2-/3-")
relations <- ibd3 %>%

  dplyr::select(seqID,seqID2,value,District)

for.checks <- ibd3 %>%
  dplyr::select(seqID,seqID2,value,District,District2)

relations <- ibd3 %>%
  dplyr::select(seqID,seqID2,value,District)

colnames(relations) <- c("from","to","F","District")

relations <- relations[order(relations$District),]
relations$District <- as.factor(relations$District)


ibd3 <- ibd2 %>%
  filter(value == 1.0) %>%
filter (hrp23_status_f == "pfhrp2-/3-")

relations <- ibd3 %>%
  dplyr::select(seqID,seqID2,value,District)

colnames(relations) <- c("target","source","F","District")
relations <- relations[order(relations$District),]

for.checks <- ibd3 %>%
  dplyr::select(seqID,seqID2,value,District,District2)
colnames(for.checks) <- c("target","source","F","District","District2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

linksdoubledeleted <- relations %>% dplyr::select(target, source, F)
#write.csv(linksdoubledeleted, "linksdoubledeleted.csv" )

seq1 <- for.checks %>% dplyr::select(target,District) 
seq2 <- for.checks %>% dplyr::select(source,District2) %>% rename(target=source,District=District2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,District)
nodes <- nodes[order(nodes$District),]
nodes <- droplevels(nodes)

#tardoubledeleted <- metadcomplete %>% 
 # filter(metadcomplete$seqID %in% nodes$target) %>%
 # dplyr::select(seqID, District)

#write.csv(tardoubledeleted, "tardoubledeleted.csv") # modify and order per pop

Fig6C_input<-read.csv("Fig6C_input.csv", header = T)

#ta2 <- ta %>%
 # dplyr::select(seqID, hrp23_status_f)


gd1 <- graph_from_data_frame(d=linksdoubledeleted, vertices=Fig6C_input, directed=F) 


#colorsforyourchoice = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")



colorsdh1 <- c(rep("#60bf37",5), rep("#953ada",12), rep("#862977",2),rep("#bba672",2),rep("#403367", 5), rep("#dd6bbb", 6), rep("#da8a6d", 4), rep("#a79cd4", 16), rep ("black", 4), rep("#d593a7", 12),rep("#895c8b", 18))

plot(gd1, vertex.label=NA,vertex.size=1, main="Only R622I Mutant IBD>=0.95 per District", vertex.color=colorsdh1) 
legend(x="topleft", legend=c("Ahferom", "At_Tsimbila", "Gulomekeda","Itang", "K_Humera", "Kule", "L_Adiabo","Metema", "Quara","Tegede", "West_Armachiho"),
      col=c( "#60bf37","#953ada","#862977","#bba672","#403367", "#dd6bbb", "#da8a6d",  "#a79cd4","black", "#d593a7","#895c8b"), cex=0.8, pch=c(19))


# Save the plot 
ggsave("Fig6C.svg", dpi=600, width=7.5, height=7)
ggsave("Fig6C.pdf", dpi=600, width=7.5, height=7)


sessionInfo()
