

#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################


####################################################
# Fig5- IBD analysis   
####################################################

getwd()
setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig5")



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
# Estimate IBD 
#######################
#######################

# Read in vcf file with only choma samples 
invisible(vcf <- read.vcfR("Only_COI1_EPHI_target_IBC_biallelicsnpmiss05sample05.recode.vcf",verbose=FALSE))

#Read in target snp positions

test <- vcf2mipanalyzer_biallelic(vcfR=vcf,verbose=FALSE)
invisible(inf <- inbreeding_mle(test, f=seq(0,1,l=100), report_progress=FALSE,ignore_het=TRUE))

# transform ibd into a data frame (var1 and var2 are the samples being compared, but referring to the index of their order in the vcf file)
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

metadcomplete_622I<-read.csv("Fig5_metadata.csv", header = T)

#check the cooridinates

metadcomplete_622I$Long <- ifelse(metadcomplete_622I$District == "Ahferom", 39.15, 
                                  ifelse(metadcomplete_622I$District == "At_Tsimbila", 38.25, metad$Lat))

metadcomplete_622I$Lat <- ifelse(metadcomplete_622I$District == "Ahferom", 14.28, 
                                 ifelse(metadcomplete_622I$District == "At_Tsimbila", 14.11, metadcomplete_622I$Lat))

###############
###############
# Fig1A - IBD distribution
###############
###############

ibd_input<-read.csv("ibd_input_coi1.csv")

mainplot2 <- ibd_input %>%
  ggplot() +
  geom_histogram(aes(x=value, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("Pairwise relatedness (IBD)") + ylab("frequency (%)") +
  theme_classic()

insetplot <- ibd_input %>%
  ggplot() +
  geom_histogram(aes(x=value, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("IBD>0.5") + ylab("frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,1.5)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
# cowplot
cowplot::ggdraw() +
  cowplot::draw_plot(mainplot2, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)


# Save plot 
ggsave("Fig5A.svg", dpi=600, width=7.5, height=7)
ggsave("Fig5A.pdf", dpi=600, width=7.5, height=7)



###############
###############
# Fig5B - R622I mutant vs wild parasite IBD comparison 

###############
###############

Fig5B<- read.csv("Fig5B_file.csv")


colors=c("#D95F02", "#1b98e0")


my_comparisons2 <- list( c("Mutant","Wild"))

ggbox2<- ggboxplot(Fig5B, x = "k13_R622I", y = "value", fill= "k13_R622I", cex=0.5)
ggbox2+
  scale_y_continuous()+
  #theme_gray()+
  scale_fill_manual(values=c("#D95F02", "#1b98e0"))+
  #theme_bw () +
  labs(x="R622I status", y="Pairwise IBD value")  +
  
  #theme(axis.text.x = element_text(angle =90))+ 
  stat_compare_means(comparisons=my_comparisons2 )+ 
  # Default method = "kruskal.test" for multiple groups # Global p-value
  ggtitle("IBD value 622I Mutant vs Wild parasites across Ethiopia")



# Save plot
ggsave("Fig5B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig5B.pdf", dpi=600, width=7.5, height=7)



###############
###############
#Fig5C - networks for clonal parasites (IBD>=0.95) per k13 R622I mutation status
###############
###############
library(dplyr)
dbs2 <- ibd_input

dbs3 <- dbs2 %>%
  filter(value==1 )
relations <- dbs3 %>%
  
  dplyr::select(seqID,seqID2,value,k13_R622I)

for.checks <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,k13_R622I,k13_R622I2)

relations <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,k13_R622I)

colnames(relations) <- c("from","to","F","k13_R622I")

relations <- relations[order(relations$k13_R622I),]
relations$k13_R622I <- as.factor(relations$k13_R622I)

dbs3 <- dbs2 %>%
  filter(value == 1 )
relations <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,k13_R622I)

colnames(relations) <- c("target","source","F","k13_R622I")
relations <- relations[order(relations$k13_R622I),]

for.checks <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,k13_R622I,k13_R622I2)
colnames(for.checks) <- c("target","source","F","k13_R622I","k13_R622I2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

linksR622I <- relations %>% dplyr::select(target, source, F)
#write.csv(links, "links.csv" )

seq1 <- for.checks %>% dplyr::select(target,k13_R622I) 
seq2 <- for.checks %>% dplyr::select(source,k13_R622I2) %>% rename(target=source,k13_R622I=k13_R622I2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,k13_R622I)
nodes <- nodes[order(nodes$k13_R622I),]
nodes <- droplevels(nodes)

#taR622I <- metadcomplete_622I %>% 
  #filter(metadcomplete_622I$seqID %in% nodes$target) %>%
  #dplyr::select(seqID, k13_R622I)


#write.csv2(taR622I, "Fig5C.csv") # modify and order per pop

Fig5C<-read.csv("Fig5C_file", header = T) # Order the lis per mutant vs wildtype

gR <- graph_from_data_frame(d=linksR622I, vertices=Fig5C, directed=F) 

colors <- c(rep("#D95F02", 22), rep("#1b98e0", 128))

# Plot 

plot(gR, vertex.label=NA,vertex.size=7, main="IBD>=0.95, n=150, Wild=128, Mutant=22", vertex.color=colors)
legend(x="topleft", legend=c("R622I Mutant","R622I Wild" ),
       col=c("#D95F02", "#1b98e0"), cex=1.2, pch=c(19))

# Save plot
ggsave("Fig5C.svg", dpi=600, width=7.5, height=7)
ggsave("Fig5C.pdf", dpi=600, width=7.5, height=7)


###############
###############
#R622I parasites per District
###############
###############

dbs2 <- ibd_input

dbs3 <- dbs2 %>%
  
  filter(grepl('Mutant', k13_R622I)) %>%
  filter (value>=0.95)

relations <- dbs3 %>%
 dplyr::select(seqID,seqID2,value,k13_R622I, District)

for.checks <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,District, District2, k13_R622I,k13_R622I2)

relations <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,District, k13_R622I)

colnames(relations) <- c("from","to","F","District")

relations <- relations[order(relations$District),]
relations$District <- as.factor(relations$District)

dbs3 <- dbs2 %>%
  
  filter(grepl('Mutant', k13_R622I)) %>%
  
  filter (value>=0.95)

relations <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,District, k13_R622I)

colnames(relations) <- c("target","source","F","District")
relations <- relations[order(relations$District),]

for.checks <- dbs3 %>%
  dplyr::select(seqID,seqID2,value,District, District2, k13_R622I,k13_R622I2)
colnames(for.checks) <- c("target","source","F","District","District2")

test <- for.checks #%>% filter(Season != Season2)
test1 <- as.data.frame(test$target)
test2 <- as.data.frame(test$source)
colnames(test2) <- "target"
colnames(test1) <- "target"
test3 <- rbind(test1,test2)
test4 <- test3 %>% distinct(target)

links <- relations %>% dplyr::select(target, source, F)
#write.csv(links, "links622I.csv" )

seq1 <- for.checks %>% dplyr::select(target,District) 
seq2 <- for.checks %>% dplyr::select(source,District2) %>% rename(target=source,District=District2)
n1 <- rbind(seq1,seq2)
nodes <- n1 %>% distinct(target,District)
nodes <- nodes[order(nodes$District),]
nodes <- droplevels(nodes)

#Fig5D <- metadcomplete_622I %>% 
  #filter(metadcomplete_622I$seqID %in% nodes$target) %>%
  #dplyr::select(seqID, District)


#write.csv(Fig5D, "Fig5D.csv") # modify and order per pop

Fig5D<-read.csv("Fig5D_file.csv", header = T)

gdh <- graph_from_data_frame(d=links, vertices=Fig5D, directed=F) 

colorsdh <- c(rep("#60bf37",1), rep("#862977",1),rep("#bba672",1),rep("#403367", 1),rep("#da8a6d", 1), rep("#a79cd4", 5), rep("#55baad", 2), rep("#d593a7", 18),rep("#895c8b"))

plot(gdh, vertex.label=NA,vertex.size=10, main="Only R622I Mutant IBD>=0.95 per District", vertex.color=colorsdh)
legend(x="topleft", legend=c("Ahferom","Gulomekeda",  "Itang", "K_Humera",  "L_Adiabo","Metema", "T_Adiado","Tegede", "West_Armachiho"),
      col=c( "#60bf37","#862977","#bba672","#403367","#da8a6d","#a79cd4","#55baad","#d593a7","#895c8b"), cex=0.5, pch=c(19))

#save plot
ggsave("Fig5D.svg", dpi=600, width=7.5, height=7)
ggsave("Fig5D.pdf", dpi=600, width=7.5, height=7)


