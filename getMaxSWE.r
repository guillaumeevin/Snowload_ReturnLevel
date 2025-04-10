# Guillaume Evin
# 09/04/2026
# guillaume.evin@inrae.fr
# maxima SWE MTMSI

library(ncdf4)
library(CFtime)

# wget -r -np -nH -m -nd -A "max-swe-NS*" -R index.html https://climatedata.umr-cnrm.fr/public/dcsc/projects/CLIMTOUR/INDICATEURS/PAR_MODELE/
fmain = "../MTMSI-max-swe/"

# vector available runs
vecruns = rbind(c("CCLM4-8-17","CLMcom-CCLM4-8-17","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP45"),
                c("CCLM4-8-17","CLMcom-CCLM4-8-17","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP85"),
                c("ALADIN53","CNRM-ALADIN53","CNRM-CM5","CNRM-CERFACS-CNRM-CM5","RCP45"),
                c("ALADIN53","CNRM-ALADIN53","CNRM-CM5","CNRM-CERFACS-CNRM-CM5","RCP85"),
                c("WRF331F","IPSL-INERIS-WRF331F","IPSL-CM5A-MR","IPSL-IPSL-CM5A-MR","RCP45"),
                c("WRF331F","IPSL-INERIS-WRF331F","IPSL-CM5A-MR","IPSL-IPSL-CM5A-MR","RCP85"),
                c("REMO2009","MPI-CSC-REMO2009","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP26"),
                c("REMO2009","MPI-CSC-REMO2009","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP45"),
                c("REMO2009","MPI-CSC-REMO2009","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP85"),
                c("RCA4","SMHI-RCA4","CNRM-CM5","CNRM-CERFACS-CNRM-CM5","RCP45"),
                c("RCA4","SMHI-RCA4","CNRM-CM5","CNRM-CERFACS-CNRM-CM5","RCP85"),
                c("RCA4","SMHI-RCA4","EC-EARTH","ICHEC-EC-EARTH","RCP26"),
                c("RCA4","SMHI-RCA4","EC-EARTH","ICHEC-EC-EARTH","RCP45"),
                c("RCA4","SMHI-RCA4","EC-EARTH","ICHEC-EC-EARTH","RCP85"),
                c("RCA4","SMHI-RCA4","IPSL-CM5A-MR","IPSL-IPSL-CM5A-MR","RCP45"),
                c("RCA4","SMHI-RCA4","IPSL-CM5A-MR","IPSL-IPSL-CM5A-MR","RCP85"),
                c("RCA4","SMHI-RCA4","HadGEM2-ES","MOHC-HadGEM2-ES","RCP45"),
                c("RCA4","SMHI-RCA4","HadGEM2-ES","MOHC-HadGEM2-ES","RCP85"),
                c("RCA4","SMHI-RCA4","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP45"),
                c("RCA4","SMHI-RCA4","MPI-ESM-LR","MPI-M-MPI-ESM-LR","RCP85"))
df.run = data.frame(rcmtag = vecruns[,1],rcm=vecruns[,2],gcmtag=vecruns[,3],
                    gcm=vecruns[,4],rcp=vecruns[,5])
saveRDS(df.run,file = "./dfRuns.rds")
nrun = nrow(df.run)

# initialize array
vyall = 1951:2100
nyears = length(vyall)
nNuts = 6584
maxSWE = array(dim=c(nrun,nNuts,nyears))

# loop over the different runs
for(irun in 1:nrun){
  dfi = df.run[irun,]
  
  # histo
  fh = paste0(fmain,"max-swe-NS-year_",dfi$rcm,
              "_",dfi$gcm,"_HISTORICAL.nc")
  nc <- nc_open(fh)
  vh = ncvar_get(nc,"max-swe-NS-year")
  vectime <- ncvar_get(nc, "time")
  tunits <- ncatt_get(nc,"time","units")
  vecTimeCf <- CFtime(tunits$value, calendar = "proleptic_gregorian", vectime) # convert time to CFtime class
  vecTimeDates <- as_timestamp(vecTimeCf, format = "timestamp")
  yh <- as.numeric(format(as.Date(vecTimeDates),"%Y"))
  nc_close(nc)
  
  # rcp
  ff = paste0(fmain,"max-swe-NS-year_",dfi$rcm,
              "_",dfi$gcm,"_",dfi$rcp,".nc")
  nc <- nc_open(ff)
  vf = ncvar_get(nc,"max-swe-NS-year")
  vectime <- ncvar_get(nc, "time")
  tunits <- ncatt_get(nc,"time","units")
  vecTimeCf <- CFtime(tunits$value, calendar = "proleptic_gregorian", vectime) # convert time to CFtime class
  vecTimeDates <- as_timestamp(vecTimeCf, format = "timestamp")
  yf <- as.numeric(format(as.Date(vecTimeDates),"%Y"))
  nc_close(nc)
  
  # filter years
  isinperiod = vyall%in%c(yh,yf)
  
  # bind histo and rcp
  maxSWE[irun,,isinperiod] = cbind(vh,vf)
}

saveRDS(maxSWE,file = "../MTMSI-max-swe/maxSWE.rds")
