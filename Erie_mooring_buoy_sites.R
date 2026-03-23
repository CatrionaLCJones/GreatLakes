library(sp)
library(sf)
library(tidyverse)
library(shapefiles)
library(viridis)
library(egg)
library(magick)

#Read in shapefile and raw data
Erie.shp <- st_read("erie_shapefile_simplified.shp")

Buoys <- read.csv("Buoy locations.csv")

#Convert the csv files to sf objects and set the CRS to match the Erie shapefile
Buoys <- st_as_sf(Buoys, coords = c("Longitude", "Latitude"), crs=4326)

Non_EPA_Buoys <- Buoys %>% filter(Source != "EPA")
  
ggplot() +
  geom_sf(data = Erie.shp) +
  geom_sf(data=Non_EPA_Buoys, aes(color=Source, shape=Location), size=8) +
  geom_sf_text(data = Non_EPA_Buoys,aes(label=Short.name))+
  xlab("Longitude") + ylab("Latitude") +
  #scale_color_viridis(discrete=TRUE)+
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) 

ggplot() +
  geom_sf(data = Erie.shp) +
  geom_sf(data=EPA, size=4, fill="green", color="black",shape=21) +
  geom_sf_text(data = EPA,aes(label=Station))+
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

#1000*450
