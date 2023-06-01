

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

setwd("C:/Users/afola/Desktop/Abefola_github/EPHI_622I_hrp23/Suplmentary_Figures/FigS1") # change this to your working directory


#Install Library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)


EPHI_metadatata<-read.csv("EPHI_622I_metadata.csv", header=T,  sep=",")


#count number of DBS samples at Region

EPHI_metadata %>% 
  group_by(Reg) %>% 
  tally()  %>% 
  arrange(n)

# Number of samples in each Region:

barplotreg <- ggplot(EPHI_metadata, aes(Reg))

barplotreg  + 
  geom_bar(aes(fill=Reg)) +
  scale_fill_brewer(palette="Dark2")+
  labs(x="Region", y="Number of samples") +
  guides(fill=guide_legend(title="Regions/Sites"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme() +
  ggtitle("Samples Distribution per Region") 

# Save plot 
ggsave("FigS1insert.svg", dpi=600, width=7.5, height=7)
ggsave("FigS1insert.pdf", dpi=600, width=7.5, height=7)

##################################################
# google map method (ggmap) and other available r package 
##################################################

setwd("C:/Users/afola/Desktop/map")

#More Q&A - https://github.com/dkahle/ggmap/issues/51

#Get the latest Install
#if(!requireNamespace("devtools")) library("devtools")
#devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)

#Set your API Key
#ggmap::register_google(key = "AIzaSyAVLxVNbGChtcW3e2gzyxOPuhnrsVbwtco")
#Use unix to get api key
#touch .Renviron
#open .Renviron

##########
# Malaria atlas -https://rspatialdata.github.io/malaria.html
#load required libraries for base maps
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(colorspace)
library(dplyr)

#2018 prevalence data (data downloaded from malaria atlas project)
prev <- read.csv("ethio_malaria_prev_data.csv")

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

ETHO<-st_read("gadm36_ETH_0.shp") # Source https://www.diva-gis.org/datadown
ETHO_regions<-st_read("gadm36_ETH_1.shp")
View(ETHO_regions)
ETHO_districts<-st_read("gadm36_ETH_2.shp") # name NAME_2 for combining
View(ETHO_districts)
#Set coordinate reference system for ETHOia
ETHO_crs<-st_crs(ETHO)
africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf", )

#water<-st_read("Africa_waterbody.shp", SHAPE_RESTORE_SHX= YES)

Malaincidence <- right_join(ETHO_regions, prev, by = "NAME_1")

#incidence

samplessites <-read.csv("ephi_lat_long.csv", header=T)
cor.reg <- c("#1B9E77", "#D95F02", "#7570B3")

pop1 = factor(samplessites$Region)

ggplot() + 
  geom_sf(data=africa, fill="gray90")+
  geom_sf(data = ETHO, fill="gray96", lwd=1.75) +
  geom_sf(data = ETHO_regions, fill="gray90") +
  geom_sf(data = Malaincidence,aes(geometry=geometry, fill=Value
  )) +
  scale_fill_distiller(type="seq", palette = "Reds", direction = 1, name="2018 Pf Incidence Rate")+
  geom_sf(data=lakes10, color="grey40", fill ="aliceblue", size= 0.8)  +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  
  #coord_sf(xlim = c(29, 41), ylim = c(-12, 0), expand = FALSE)+xlab("Longitude")+ ylab("Latitude")+
  # coord_sf(xlim = c(22.3, 33.3), ylim = c(-8,-18), expand = TRUE) +
  coord_sf(xlim = c(33, 48), ylim = c(2,15), expand = F) +
  ggtitle(expression(paste("P.falciparum Infection Prevalence	per 100 Children")))+
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  xlab("") +
  
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(data=samplessites, position=position_jitter(width=0.0, height=0.0), size = 4, pch=19, col=cor.reg[as.integer(pop1)], 
              # aes(x=lat, y=long),alpha = 6/10)+ #, color=variety, sha5e=factor(status)
              # aes(x=lat, y=long, alpha = 6/10, color=variety, shape=factor(status))) + ### if you want to color and shape
              aes(x=Longitude, y=Latitude, col= "Region" )) 

#prev Infection Prevalence	per 100 Children

ggplot() + 
  geom_sf(data=africa, fill="gray90")+
  geom_sf(data = ETHO, fill="gray96", lwd=1.75) +
  geom_sf(data = ETHO_regions, fill="gray90") +
  geom_sf(data = Malaincidence,aes(geometry=geometry, fill=Value
  )) +
  scale_fill_distiller(type="seq", palette = "PuRd", direction = 1, name="Pf-prevalence/100 Children")+
  geom_sf(data=lakes10, color="grey40", fill ="aliceblue", size= 0.8)  +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  
  #coord_sf(xlim = c(29, 41), ylim = c(-12, 0), expand = FALSE)+xlab("Longitude")+ ylab("Latitude")+
  # coord_sf(xlim = c(22.3, 33.3), ylim = c(-8,-18), expand = TRUE) +
  coord_sf(xlim = c(33, 48), ylim = c(2,15), expand = F) +
  ggtitle(expression(paste("P.falciparum Infection Prevalence	per 100 Children")))+
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  xlab("") +
  
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(data=samplessites, position=position_jitter(width=0.0, height=0.0), size = 4, pch=19, col=cor.reg[as.integer(pop1)], 
              # aes(x=lat, y=long),alpha = 6/10)+ #, color=variety, sha5e=factor(status)
              # aes(x=lat, y=long, alpha = 6/10, color=variety, shape=factor(status))) + ### if you want to color and shape
              aes(x=Longitude, y=Latitude, col= "Region" )) 

# Save plot 
ggsave("P.falciparum prevalence per region.svg", dpi=600, width=7.5, height=7)
ggsave("P.falciparum prevalence Rate per region.pdf", dpi=600, width=7.5, height=7)


