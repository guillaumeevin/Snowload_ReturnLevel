# Guillaume Evin
# 27/08/2024
# guillaume.evin@inrae.fr
# 
# produce a map showing what is the best distribution

# load libraries
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(scales)#squish
library(tidyterra)
library(ggmap)

#===============================================================================
# Shape nuts file
#===============================================================================

# read shape files
shp_6584nuts <- st_read(file.path("../NUTS/NUTS3wtAlt.shp")) # 6584 NUTS
n.nuts.alt = nrow(shp_6584nuts)
shp_nuts_raw <- st_read(file.path("../NUTS/nut_elev.shp")) # 1522 NUTS
shp_nuts = shp_nuts_raw[shp_nuts_raw$NUTS_ID%in%shp_6584nuts$NUTS_ID,] # 1515 NUTS
n.nuts = nrow(shp_nuts)

# country borders #
worldmap <- ne_countries(scale = 'medium',type = 'map_units',returnclass = 'sf')
worldmap3857 <- st_transform(worldmap, 3857) # we convert crs to pseudo mercator

# statistics for the 1517 nuts
nuts.stats = read.csv("../NUTS/uerra5_alt_ideal_mean_nuts.csv")

# world shapes
world <- ne_countries(scale = "medium", returnclass = "sf")

# build data.frame with unique nuts/altitude combination, based on the nuts with
# the altitude which are the closer to the real altitude
id = shp_nuts$NUTS_ID # 1515
dfNuts = data.frame(id=id,alt=vector(length=n.nuts),i6584=vector(length=n.nuts))
for(inuts in 1:n.nuts){
  id.i = dfNuts$id[inuts]
  
  # real altitude
  if(id.i=="NO060"){
    i.stat = which(nuts.stats$nuts_id=="NO061")
    alt.stat = nuts.stats$alt_mean[i.stat]
  }else{
    i.stat = which(nuts.stats$nuts_id==id.i)
    alt.stat = nuts.stats$alt_mean[i.stat]
  }
  
  # find corresponding nut
  indices6584 = which(shp_6584nuts$NUTS_ID==id.i)
  alts6584 = shp_6584nuts$alt[indices6584]
  isel = which.min(abs(alts6584-alt.stat))
  
  dfNuts$alt[inuts] = alts6584[isel]
  dfNuts$i6584[inuts] = indices6584[isel]
}

#===============================================================================
# read distribution codes
#===============================================================================
n.couple = 9 # number of GCM/RCM couples
best.dist = vector(length = n.nuts.alt)
for(inut in 1:n.nuts.alt){
  f = paste0("../ReturnLevelOutput/l",inut,".rds")
  if(file.exists(f)){
    # list of models outputs for the 9 GCM?RCM couples
    ll = readRDS(file = f)
    
    # vector of the best models
    vecbest = ll$ibest
    
    # count best models
    countbest = table(vecbest)
    ibest = which(countbest==max(countbest))
    
    if(length(ibest)==1){
      bestmodels = names(countbest)
      best.dist[inut] = bestmodels[ibest]
    }else{
      best.dist[inut] = "Mixed"
    }
  }else{
    best.dist[inut] = "no fit"
  }
}
best.dist.map = best.dist[dfNuts$i6584]

plt =  ggplot() +
  geom_sf(data=worldmap)+
  geom_sf(data = shp_nuts, aes(fill = best.dist.map), color = NA) +
  guides(fill=guide_legend(title="",reverse = TRUE)) +
  coord_sf(xlim = c(-30, 50),
           ylim = c(34, 71))+
  theme_light()+
  theme(axis.text = element_text(size=6))+
  scale_fill_manual(breaks=c("0","Mixed","GEV","GUM","GG","GA","pIG","pGA","pEXP"), 
                      labels=c(">90% of zeros","No dominant model","GEV","GUM","GG","GA","pIG","pGA","pEXP"),
                    values = c("gray20","darkkhaki","blue","cyan","magenta4","magenta1","orange","red","yellow"))
ggsave(plot = plt, filename = "Fig6_mapbestmodel_raw.pdf", path = "../Figures/",
       device = "pdf",units="cm",height=12,width=18)