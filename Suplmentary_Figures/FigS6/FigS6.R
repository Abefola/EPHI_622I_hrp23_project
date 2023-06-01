
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

#################
#################
# FigS6A - DHPS-DHFR UpSet plot quantiple or setplex mutants
#################
#################

getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS6") # change this to your working directory


#Install Library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(UpSetR)

#Load the data 

FigS6A <- read.csv("FigS6A_file.csv", header=TRUE, sep="," )


#Columns contain
# Sample Id, mutant names, and Wild genotype (if it contains wild alleles for all loci )

### Define colors 

bar_cols1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
bar_cols <- c("#0072B2","darkgreen","darkcyan", "grey", "black", "red","tomato4","orange")

#bar_cols <- c( "#E69F00", "#CC79A7", "#009E73","#F0E442", "#0072B2", "#D55E00",  "#56B4E9", "black")

upset(FigS6A, nsets = 9, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE), keep.order = FALSE, main.bar.color = bar_cols1)

# Save plot 
ggsave("FigS6A.svg", dpi=600, width=7.5, height=7)
ggsave("FigS6A.pdf", dpi=600, width=7.5, height=7)



#################
#################
# FigS6B special prevalence of 581G
#################
#################

#laod requered libraries
library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(scatterpie)
library(scales)
library (ggforce)
library (ggmap)
library (ggplot2)

#https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/

#load base map layers
admin8 <- ne_download(scale="large", type = "admin_1_states_provinces_lines",
                      category = "cultural", returnclass = "sf")
rivers8 <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                       category = 'physical', returnclass = "sf")
lakes8 <- ne_download(scale = "large", type = 'lakes',
                      category = 'physical', returnclass = "sf")
sov10 <- ne_download(scale="large", type = "sovereignty",
                     category = "cultural", returnclass = "sf")
admin10<- ne_download(scale="large", type = "populated_groups",
                      category = "cultural", returnclass = "sf")

# Load the the data (dhps A581G_per_district) 

FigS6B <- read.csv("FigS6B.csv", header=T)

#scale pies by sample size
FigS6B$radius <- rescale(FigS6B$total, to = c(2, 5))


#Plot

ggplot()+
  geom_sf(data=sov10, color='black', size=1.5,
          fill = ifelse(sov10$ADMIN == "Ethiopia", '#BBFFFF', 'grey90')) +
  #geom_sf(data=rivers8, color="cyan4", size=0.5, alpha=0.5) +
  geom_sf(data=lakes8, color="grey40", fill ="aliceblue", size= 0.8) +
  geom_sf(data=admin8, color="grey", size= 1.5) +
  geom_sf(data=admin10, color="grey80", size= 0.1) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(which_north = "true",  height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         pad_x = unit(0.85, "cm"),
                     pad_y = unit(0.6, "cm"))+
  
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in")) +
  #style = north_arrow_fancy_orienteering)+
  annotate("text", x = 34.5, y = 13, label = "Sudan",
           color="grey30", size=5 , fontface="italic") +
  annotate("text", x = 32.5, y = 9.5, label = "South Sudan", angle= 90,
           color="grey30", size=5 , fontface="italic") +
  annotate("text", x = 37.5, y = 15, label = "Eritrea",
           color="grey30", size=5 , fontface="italic") +
  geom_scatterpie(aes(x=long, y=lat, r = (radius/18)),
                  data = FigS6B,
                  cols = c("A581G","WT"),
                  color = "grey20",
                  legend_name = "Pfdhps_A581G",
                  alpha=1)+
  scale_fill_manual(values = c("red","white")) +
  coord_sf(xlim = c(31.6, 40), ylim = c(7.5,15), expand = T) +
  theme_void()+
  theme(panel.background = element_rect(fill = 'aliceblue', colour = "grey80"))

# Save plot 
ggsave("FigS6B.svg", dpi=600, width=7.5, height=7)
ggsave("FigS6B.pdf", dpi=600, width=7.5, height=7)

