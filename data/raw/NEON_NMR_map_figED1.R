# make map of sampling locations (Figure S1)

#remove old data
rm(list=ls())
library(ggmap)
library(ggplot2)
library(ggrepel)

data<-read.csv("Hall_Ye_2020_NEON_NMR_20200716.csv",header=TRUE)


# make a map of US + PR
map.usa<-get_stamenmap(bbox=c(left= -127,bottom = 16, right=-64, top =50),
                        zoom=4,
                  maptype="terrain-background")

map.ak<-get_stamenmap(bbox=c(left= -173,bottom = 55, right=-140, top =71),
                       zoom=4,
                       maptype="terrain-background")

colors<-viridis(6)

set.seed(1)

quartz(,5,3.75)
ggmap(map.usa)+
  # theme(axis.line=element_blank(), # get rid of axis notation
  #       axis.text.x=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks=element_blank(),
  #       axis.title.x=element_blank(),
  #       axis.title.y=element_blank())+
  geom_point(data = data, # plot SOC constituents
             aes(x = Longitude_degrees, 
                 y = Latitude_degrees),
             alpha=0.7,
             colour="darkblue",
             size=1)+
  geom_text_repel(data = data, # plot SOC constituents
             aes(x = Longitude_degrees, 
                 y = Latitude_degrees,
                 label=siteID),
             alpha=0.7,
             colour="darkblue",
             size=3)
quartz(,1.5,1.5)
ggmap(map.ak)+
  theme(
        
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  geom_point(data = data, # plot SOC constituents
             aes(x = Longitude_degrees, 
                 y = Latitude_degrees),
             alpha=0.7,
             colour="darkblue",
             size=1)+
  geom_text_repel(data = data, # plot SOC constituents
                  aes(x = Longitude_degrees, 
                      y = Latitude_degrees,
                      label=siteID),
                  alpha=0.7,
                  colour="darkblue",
                  size=3)
