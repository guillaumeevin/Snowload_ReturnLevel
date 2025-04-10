# ---
# title: "Mapping elevation"
# output: html_document
# date: "2024-07-03"
# ---

## Library
library(tidyverse)
library(rnaturalearth)
library(tidyr)
library(sf)
library(ggmap)
library(ggplot2)
library(plyr)
library(dplyr)
library(sf)
library(RColorBrewer)
library(Polychrome)
library(ggrepel)


## Mapping
#
#Ajusted elevation = closest elevation band (from MTMSI) to ideal mean elevation of the NUTS

# LOADING DATA #
shp_alt_ajust3857 <- st_read(file.path('.',"shp_alt_ajuste_3857.shp")) # shapefile of ajusted elevation per nuts

# COLOUR PALETTE #
absrl_max=max(shp_alt_ajust3857$alt_ajust, na.rm=TRUE)
absrl_min=min(shp_alt_ajust3857$alt_ajust, na.rm=TRUE)
limite=c(absrl_min, absrl_max)
palette_alti = c("#08519C", "#427FE2", "#B2E2E2", "#41AB5D", "#C2E699",
                 "khaki", "#FEC44F", "#FE9929", "#CC4C02", "#993404", "#662506")

# MAPPING #

# country borders #
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
worldmap3857 <- st_transform(worldmap, 3857) # we convert crs to pseudo mercator

# mapping with ggplot #
nut1 = shp_alt_ajust3857[shp_alt_ajust3857$nuts_id=="AT333",]
nut2 = shp_alt_ajust3857[shp_alt_ajust3857$nuts_id=="FRG01",]

plt = ggplot() +
  geom_sf(data=worldmap)+
  geom_sf(data=shp_alt_ajust3857,aes(fill = alt_ajust), colour=NA)+
  scale_fill_gradientn(colours=palette_alti, limits = limite, na.value = NA, 
                       guide = guide_colorbar(barwidth = 0.5, barheight = 10, nbin = 11, title="Elevation (m)"))+
  geom_sf(data = nut1,fill = "cyan",colour="black",linewidth=0.4)+
  geom_sf(data = nut2,fill = "cyan",colour="black",linewidth=0.4)+
  coord_sf(xlim = c(-30, 50),
           ylim = c(34, 71))+
  theme_light()+
  theme(axis.text = element_text(size=6))


ggsave(filename = "../Figures/Fig1_mapElevNUTS.pdf",
       plot=plt,device = "pdf",units="cm",height=12,width=18)