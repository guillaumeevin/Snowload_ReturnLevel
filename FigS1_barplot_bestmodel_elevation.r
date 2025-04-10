# Guillaume Evin
# 10/04/2025
# guillaume.evin@inrae.fr
# 
# barplot of the number of times models are selected against elevation

# load libraries
library(sf)
library(colorspace)
library(ggplot2)
library(scales)
library(ggh4x)

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
# read distribution codes
#===============================================================================
n.couple = 9 # number of GCM/RCM couples
best.dist = matrix(nrow = n.nuts.alt, ncol=n.couple)
for(inut in 1:n.nuts.alt){
  f = paste0("../ReturnLevelOutput/l",inut,".rds")
  if(file.exists(f)){
    # list of models outputs for the 9 GCM?RCM couples
    ll = readRDS(file = f)
    
    # vector of the best models
    best.dist[inut,] = ll$ibest
  }
}

# convert as a data.frame
df.raw = data.frame(dist=factor(as.vector(best.dist[dfNuts$i6584,]),
                                levels = c("pEXP","pGA","pIG","GA","GG","GUM","GEV","0")),
                alt=as.vector(replicate(n = n.couple,dfNuts$alt)))
df = df.raw[df.raw$dist!="0",]

design <- "
 ABC#
 DEFG
"

plt = ggplot(df, aes(x=alt)) + 
  geom_bar(stat = "count", position = "dodge2")  +
  facet_manual(~dist,design) + scale_y_continuous(trans = "log10") +
  labs(x = "Elevation (m)",y="Number of model selections") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill=NA),
    strip.text = element_text(face="bold")
  )

ggsave(filename = paste0("FigS1_barplot_bestmodel_elevation.pdf"),
       plot =  plt, path = "../Figures/",
       device = "pdf",units="cm",height=10,width=18)