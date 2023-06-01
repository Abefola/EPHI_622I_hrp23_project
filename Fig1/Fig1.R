
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

#################
#################
# Part I - not weighted 622I prevalence 
#################
#################


getwd() # To check the current dir

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Fig1") # change this to your working directory

#################
# load all required libraries 
#################

#devtools::install_github("ropensci/rnaturalearth")
#remove.packages (rnaturalearthdata, lib)

library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(scatterpie)
library(scales)
library (ggforce)
library (ggmap)
library(dplyr)
library (ggplot2)
library(ggpubr)

###################
#Drug resistance spatial spread at District level - Fig 1A
###################

#readings
#https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/

############
# load base map layers 
############

admin8 <- ne_download(scale="large", type = "admin_1_states_provinces_lines",
                      category = "cultural", returnclass = "sf")
rivers8 <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                       category = 'physical', returnclass = "sf")
lakes8 <- ne_download(scale = "large", type = 'lakes',
                      category = 'physical', returnclass = "sf")
sov10 <- ne_download(scale="large", type = "sovereignty",
                     category = "cultural", returnclass = "sf")
admin10<- ne_download(scale="large", type = "populated_places",
                     category = "cultural", returnclass = "sf")
#############
#scale pies by sample size
#############

Fig1A <- read.csv("Fig1A_file.csv", header=T) # 

Fig1A$radius <- rescale(Fig1A$total, to = c(3, 5))

#############
#plot prevalence pie chat on Ethiopia map
#############
ggplot()+
geom_sf(data=sov10, color='black', size=1.5, fill = ifelse(sov10$ADMIN == "Ethiopia", 
                                                           '#BBFFFF', 'grey')) +
 # geom_sf(data=rivers8, color="cyan4", size=0.5, alpha=0.5) +
  geom_sf(data=lakes8, color="grey40", fill ="aliceblue", size= 0.8) +
  geom_sf(data=admin8, color="grey80", size= 1.5) +
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
                  data = Fig1A,
                  cols = c("Mutant","Wild"),
                  color = "grey20",
                  legend_name = "K13_622I",
                  alpha=0.7)+
  scale_fill_manual(values = c("#D95F02", "#1b98e0")) +
  #coord_sf(crs = 20138) +
  coord_sf(xlim = c(31.6, 40), ylim = c(7.5,15), expand = T) +
  theme_void()
  #theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97'))

# Save plot 
ggsave("Fig1.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1.pdf", dpi=600, width=7.5, height=7)



###################
# 622I bar plot at Regional- Fig 1A-insert bar plot 
###################

# plot 
Fig1A_insert<-read.csv('Fig1A_insert_bar_file.csv')

barplotreg<-ggplot(data=Fig1A_insert, aes(x=Region, y=Prevalence)) 

barplotreg  + 
  geom_bar(stat="identity", fill= "#D95F02")+
  #scale_fill_brewer(palette="Dark2")+
  labs(x="Region", y="Prevalence (%)") +
  theme(text = element_text(size = 20)) 
  #coord_flip() +
  #guides(fill=guide_legend(title="Regions/Sites"))+
  #theme(axis.text.x = element_text(angle = 90))+
  #theme_bw()

#save plot
# Save plot 
ggsave("Fig1Ainsert1.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1Ainsert1.1.pdf", dpi=600, width=7.5, height=7)



###############
###############
#Fig1AinsertB - Africa and Ethiopia with different color
###############
###############

#load required libraries for base maps
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(colorspace)
library(dplyr)

#get admin shape file from naturalearth
admin10 <- ne_download(scale="large", type = "admin_1_states_provinces",
                       category = "cultural", returnclass = "sf")
rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes', 
                       category = 'physical', returnclass = "sf")
oceans10 <- ne_download(scale = "medium", type = "coastline",
                        category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")

#get just the ETHOia  data # Source https://www.diva-gis.org/datadown

ETHO<-st_read("gadm36_ETH_0.shp") 
ETHO_regions<-st_read("gadm36_ETH_1.shp")
View(ETHO_regions)
ETHO_districts<-st_read("gadm36_ETH_2.shp") # name NAME_2 for combining
View(ETHO_districts)
#Set coordinate reference system for ETHOia
ETHO_crs<-st_crs(ETHO)
africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf", )
# Africa vs Eth Map

ggplot() + 
  geom_sf(data=africa, fill="gray90")+
  geom_sf(data = ETHO, fill="#BBFFFF", lwd=1.75) +
  geom_sf(data = ETHO_regions, fill="#BBFFFF") +
  theme_bw()

ggsave("Fig1AinsertB.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1AinsertB.pdf", dpi=600, width=7.5, height=7)

###################
# Fig1B Kelch mutations per WHO category 
###################

Fig1B <- read.csv("Fig1B_file.csv", header=TRUE, sep="," )

# plot
cbp2 <- c("#D95F02","orange","white") # my color choice 

kelch13 <- ggplot(Fig1B, aes(x=reorder(Mutation, AA_coden_position), y=prevalence, fill = WHO_Category
))
kelch13  + 
  geom_bar(stat = "identity")+
  labs(x="Mutation name", y="Prevalence (%)") +
  scale_fill_discrete(name="WHO Category") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
  scale_fill_manual(values=cbp2)
ggtitle("NS mutations across Kelch13 gene in Ethiopia")

# save plot

# Save plot 
ggsave("Fig1B.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1B.pdf", dpi=600, width=7.5, height=7)




###################
# Fig1C Weighted prevalence of key mutations
###################

#plot 
Fig1C<- read.csv("Fig1C_file.csv") 

kelch13 <- ggplot(keygenes, aes(x=reorder(Mutation_Name
, Gene),y=weighted_prevalence, fill = Gene))
kelch13  + 
  geom_bar(stat = "identity")+
  labs(x="Mutation name", y="Prevalence (%)") +
  scale_fill_discrete(name="Gene name") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
  #scale_fill_manual(values=cbp1)+
  ggtitle("Key mutations across three regions in Ethiopia")

#save plot

ggsave("Fig1C.svg", dpi=600, width=7.5, height=7)
ggsave("Fig1D.pdf", dpi=600, width=7.5, height=7)

