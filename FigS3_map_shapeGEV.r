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
library(latex2exp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Included here as it is not exported from ggplot2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
binned_pal <- function (palette) {
  function(x) {
    palette(length(x))
  }
}


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
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
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
# read return levels
#===============================================================================
n.couple = 9 # number of GCM/RCM couples
xi.mean = rep(NA,n.nuts.alt)
for(inut in 1:n.nuts.alt){
  f = paste0("../ReturnLevelOutput/l",inut,".rds")
  ll = readRDS(file = f)
  if(!(any(is.na(ll$vec.shape.GEV)))){
    if(any(ll$ibest=="GEV")){
      is.GEV = ll$ibest=="GEV"
      xi.vec = ll$vec.shape.GEV[is.GEV]
      xi.mean[inut] = mean(xi.vec)
    }
  }
}

#===============================================================================
# colorscale
#===============================================================================

# selected nuts/alt
shp_nuts$xi = xi.mean[dfNuts$i6584]

mybreaks <- c(-0.3,-0.2,-0.1,0,0.1,0.2)
mylabels = mybreaks

vec.col = c('#01665e','#5ab4ac','#c7eae5','#f6e8c3','#d8b365')

colorscale = binned_pal(scales::manual_pal(vec.col))

plt =  ggplot() +
  geom_sf(data=worldmap)+
  geom_sf(data = shp_nuts, aes(fill = xi), color = NA) +
  coord_sf(xlim = c(-30, 50), ylim = c(34, 71))+
  theme_light()+
  theme(axis.text = element_text(size=6))+
  guides(fill=guide_colorbar(barwidth = 1.5, barheight = 15,
                             labbel.theme = element_text(size = 10, face = c("bold"),color=c("black")),
                             title.theme=element_text(size = 10, face = "bold"),title.position = "top"))+
  binned_scale(aesthetics = "fill",name=TeX("$\\xi$"),palette = colorscale,
               guide="coloursteps",limits=range(mybreaks),breaks=mybreaks,
               show.limits = T,oob=squish,labels= mylabels)
ggsave(filename = paste0("FigS3_map_xi.png"), path = "../Figures/",
       device = "png",units="cm",height=12,width=18,dpi=300)