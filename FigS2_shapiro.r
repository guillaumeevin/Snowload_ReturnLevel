# Guillaume Evin
# 28/08/2024
# guillaume.evin@inrae.fr

library(sf)

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

i1515NUTS = dfNuts$i6584

#===============================================================================
# read p-values of the shapiro test
#===============================================================================
n.couple = 9 # number of GCM/RCM couples
mat.pval = mat.bestmodel = matrix(nrow = n.nuts,ncol = n.couple)
for(i in 1:n.nuts){
  inut = i1515NUTS[i]
  
  f = paste0("../ReturnLevelOutput/l",inut,".rds")
  
  # list of models outputs for the 9 GCM RCM couples
  ll = readRDS(file = f)
  
  # vector of the best models
  mat.pval[i,] = ll$pval
  mat.bestmodel[i,] = ll$ibest
}

labsub = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)")
pdf(file = "../Figures/FigS2_shapiro.pdf",width = 8,height = 6)
par(mfrow=c(3,3),mar=c(3.8,3.8,3,1))
isub = 0
for(labmodel in c("pEXP","pGA","pIG","GA","GG","GUM","GEV")){
  isub = isub+1
  ismodel = mat.bestmodel==labmodel
  h=hist(mat.pval[ismodel],breaks=20,
       main = paste0(labsub[isub]," ",labmodel,": n=",sum(ismodel)),
       probability = T,xlab="",ylab="",cex.lab=1.3,ylim=c(0,2.5))
  mtext(side = 1,text = "p-value",line=2.5,cex=0.8)
  mtext(side = 2,text = "Density",line=2.5,cex=0.8)
  abline(h=1,lwd=2,lty=2)
  rect(xleft = h$mids[1]-0.025,xright = h$mids[1]+0.025,ybottom = 0,ytop = h$density[1],col="red")
}
dev.off()


#===============================================================================
# read data to illustrate some bad fits
#===============================================================================

# period for the models (global temperatures available until 2099)
vy = 1951:2099

maxSWE = readRDS("../MTMSI-max-swe/MMEmaxSWE.rds")
dfRuns = readRDS("./dfRuns.rds")
source("D:/BACKUP/SnowLoad/lib_fitSWE.r")
n.scen = dim(maxSWE)[1]
n.nutsalt = dim(maxSWE)[2]
sl.years = 1951:2100 # checked
zz.sl = sl.years%in%vy


# load global temperature: vecYearsGT,HadCRUT5.spline,GCM.spline,globalTempRef,scen_GCM_RCP,scenAvail,gtas
load("D:/BACKUP/SnowLoad/GlobalTempSmooth.RData")
zz.gt = vecYearsGT%in%vy

# format a matrix 200 * 20 to get global temperatures from 1900 to 2099 (200 years)
# for each of the 20 scenarios
matGT = matrix(nrow=length(vy),ncol=n.scen)
for(is in 1:n.scen){
  iGCMtag = which(scenAvail$GCM==dfRuns$gcm[is])[1]
  GCMtag = scenAvail$GCMtag[iGCMtag]
  comb = paste0(GCMtag,"_",tolower(dfRuns$rcp[is]))
  wl = globalTempRef[[comb]]
  matGT[,is] = wl[zz.gt]
}

# get couple GCM/RCM
gcm.rcm = paste0(dfRuns$gcmtag," / ",dfRuns$rcmtag)
couple.gcm.rcm = unique(gcm.rcm)
n.couple = length(couple.gcm.rcm)

# threshold for the SWE (mm)
threshold = 1

#_______________________________________________________________________________
#  illustrate some bad fits
#_______________________________________________________________________________

iLowPval = which(mat.pval<10^-4,arr.ind = T)

nBadFits = nrow(iLowPval)

par(mfrow=c(nBadFits,2))

for(ibad in 1:nBadFits){
  inut = i1515NUTS[iLowPval[ibad,1]]
  i.couple = iLowPval[ibad,2]
  
  maxSWEnut = readRDS(file = paste0("../MMENUTS/maxSWE",inut,".rds")) # 20 scenarios * 149 years
  pred.gt = seq(from=0,to=6,by=0.1)
  
  # get ony the scenario RCP85 for that GCM/RCM
  is = which(gcm.rcm==couple.gcm.rcm[i.couple]&dfRuns$rcp=="RCP85")
  gt = as.vector(matGT[,is])
  maxSWE.couple = maxSWEnut[is,]
  
  # vectorize
  maxswe = as.vector(t(maxSWE.couple))
  
  # fit and find the best model
  df = data.frame(gt=gt,maxswe=maxswe)
  df.zz = get.df.nona(df,threshold)
  lmodel = fit.SWE(df.zz)
  i.best = i.best.model(lmodel)
  
  selmod = lmodel[[i.best]]
  res = residuals(selmod)
  shapiro.out = shapiro.test(res)
  pval = shapiro.out$p.value
  
   
  plot(df$gt,df$maxswe,pch=20,
       xlab="Global warming levels [?C]",ylab="Return levels of SWE [mm]",
       main=paste0(shp_6584nuts$NUTS_ID[inut]," / ",shp_6584nuts$alt[inut],
                   "m / ",selmod$label))
  lines(pred.gt,get.predicted.rl(selmod,pred.gt=pred.gt,proba=0.5,df.zz),col="red",lwd=2)
  lines(pred.gt,get.predicted.rl(selmod,pred.gt=pred.gt,proba=0.75,df.zz),col="blue",lwd=2)
  lines(pred.gt,get.predicted.rl(selmod,pred.gt=pred.gt,proba=0.9,df.zz),col="cyan",lwd=2)
  lines(pred.gt,get.predicted.rl(selmod,pred.gt=pred.gt,proba=0.98,df.zz),col="green",lwd=2)
  
  plot(density(res),xlim=c(-4,4),lwd=2,
       xlab="Normalized residuals",main="")
  x = seq(from=-4,to=4,by=0.01)
  lines(x,dnorm(x),col="red",lwd=2)
  
  
}
