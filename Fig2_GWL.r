# ---
#   title: "global warming levels"
# author: "Guillaume Evin"
# date: "06/08/2024"
# output: html_document
# ---


# ensemble of projections for the MTMSI dataset
dfRuns = readRDS("dfRuns.rds")

# load global temperature: vecYearsGT,HadCRUT5.spline,GCM.spline,globalTempRef,scen_GCM_RCP,scenAvail,gtas
load("../GlobalTAS/GlobalTempSmooth.RData")

#===============================================================================
## Global temperatures as a function of time
#===============================================================================

# colors GCM
library(ggthemes)
col8 = colorblind_pal()(8)

# short names
vecGCM = unique(scenAvail$GCMtag[scenAvail$GCM%in%dfRuns$gcm]) # 5 GCMs
mycol.GCM = list("CNRM-CM5"=col8[7],"EC-EARTH"=col8[8],"HadGEM2-ES"=col8[4],
                 "IPSL-CM5A-MR"=col8[2],"MPI-ESM-LR"=col8[6])

# the RCP8.5 scenario only
rcp = "rcp85"

pdf("../Figures/Fig2_GWL.pdf",width=9,height=6)
par(mar=c(4.5,5.5,0.5,0.5))

plot(-1,-1,xlim=range(vecYearsGT),ylim=c(-0.1,6),
     xlab="Years",ylab="", cex.axis=1.3,cex.lab=1.5)

# GCM
for(iGCM in 1:length(vecGCM)){
  gcm = vecGCM[iGCM]
  comb = paste0(gcm,"_",rcp)
  spl = GCM.spline[[comb]]
  lines(vecYearsGT,spl$yin-spl$y[1],lty=2,lwd=1,col=mycol.GCM[[gcm]])
  lines(vecYearsGT,globalTempRef[[comb]],lty=1,lwd=3,col=mycol.GCM[[gcm]])
}

# HadCRUT
lines(HadCRUT5.spline$x,HadCRUT5.spline$yin-HadCRUT5.spline$y[1],col="black",lty=2,lwd=1)
lines(HadCRUT5.spline$x,HadCRUT5.spline$y-HadCRUT5.spline$y[1],lty=1,lwd=3)

grid()

mtext(side=2,text = "Anomalies of global mean temperature", line=4,cex=1.5)
mtext(side=2,text = "with respect to pre-industrial levels (°C)",cex=1.5,line=2.5)

legend("topleft",bty="n",legend=c("HadCRUT5",names(mycol.GCM)),seg.len=4,cex=1.2,
       col=c("black",unlist(mycol.GCM)),lty=1,lwd=3)

legend("left",bty="n",legend=c("Smoothed global mean","Raw global mean"),
       seg.len=4,cex=1.2, col="black",lty=c(1,2),lwd=3)

dev.off()