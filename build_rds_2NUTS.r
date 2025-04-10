# ---
#   title: "global warming levels"
# author: "Guillaume Evin"
# date: "09/04/2026"
# output: html_document
# ---

library(sf)

# ensemble of projections for the MTMSI dataset
dfRuns = readRDS("dfRuns.rds")
dfRunsRCP85 = dfRuns[dfRuns$rcp=="RCP85",]

# load global temperature: vecYearsGT,HadCRUT5.spline,GCM.spline,globalTempRef,scen_GCM_RCP,scenAvail,gtas
load("../GlobalTAS/GlobalTempSmooth.RData")

# maxSWE
maxSWE = readRDS("../MTMSI/maxSWE.rds")
shp_6584nuts <- st_read(file.path(".","NUTS3wtAlt.shp")) # 6584 NUTS

# period for the models (global temperatures available until 2099)
vy = 1951:2099
zz.gt = vecYearsGT%in%vy
sl.years = 1951:2100 # checked
zz.sl = sl.years%in%vy

# short names
combgcmrcm = paste0(dfRunsRCP85$gcmtag," / ",dfRunsRCP85$rcmtag)



#_______________________
# plot 1
#_______________________
inut=which(shp_6584nuts$NUTS_ID=="FRG01") # 0m/Loire-Atlantique
maxSWEnut = maxSWE[,inut,zz.sl]
l = list()

for(i in 1:length(combgcmrcm)){
  comb = combgcmrcm[i]
  
  scen = dfRunsRCP85[i,]
  is = as.numeric(rownames(scen))
  
  wlAllYears = globalTempRef[[paste0(scen$gcmtag,"_rcp85")]]
  wl = wlAllYears[zz.gt]
  
  mx = maxSWEnut[is,]
  
  l[[comb]] = data.frame(year=vy,gwl=wl,maxSWE=mx)
}
saveRDS(l,file="./NUTSLoireAtlantique0m.rds")


#_______________________
# plot 2
#_______________________

inut=which(shp_6584nuts$NUTS_ID=="AT333"&shp_6584nuts$alt==2000) # 2000m/Osttirol/Austria
maxSWEnut = maxSWE[,inut,zz.sl]

for(i in 1:length(combgcmrcm)){
  comb = combgcmrcm[i]
  
  scen = dfRunsRCP85[i,]
  is = as.numeric(rownames(scen))
  
  wlAllYears = globalTempRef[[paste0(scen$gcmtag,"_rcp85")]]
  wl = wlAllYears[zz.gt]
  
  mx = maxSWEnut[is,]
  
  l[[comb]] = data.frame(year=vy,gwl=wl,maxSWE=mx)
}
saveRDS(l,file="./NUTSOsttirol2000m.rds")