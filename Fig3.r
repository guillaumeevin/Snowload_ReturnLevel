# ---
# title: "Figure 3: SWE maxima as a function of the global warming level, for 
# the different simulation chains and the two illustrative NUTS-3 areas. (a) 
# Loire-Atlantique, France (b) Osttirol, Austria."
# author: "Guillaume Evin"
# date: "02/12/2024"
# ---

# This script produces the Figure 3 in the article entitled "Estimating changes 
# in extreme snow load conditionally to global warming levels" by Evin. G.,
# Le Roux, E., Kamir, E. and Morin, S., submitted to Cold Regions Science and 
# Technology.

# colors
library(ggthemes)
col8 = colorblind_pal()(8)
mycol.chain = list("MPI-ESM-LR / CCLM4-8-17"="skyblue2",
                 "MPI-ESM-LR / REMO2009" = "blue1",
                 "MPI-ESM-LR / RCA4" = "darkblue",
                 "CNRM-CM5 / ALADIN53"="#CC8C3C",
                 "CNRM-CM5 / RCA4"="chocolate4",
                 "IPSL-CM5A-MR / WRF331F"="orange",
                 "IPSL-CM5A-MR / RCA4"="orangered",
                 "EC-EARTH / RCA4"="violet",
                 "HadGEM2-ES / RCA4"= "olivedrab")
vec.col = unlist(mycol.chain)
vec.chain = names(mycol.chain)

###################
# start plot
###################
pdf("../Figures/Fig3_SWEmax_2NUTS.pdf",width=9,height=12)
par(mar=c(4.5,5.5,4.5,0.5),mfrow=c(2,1))

#_______________________
# (a) Loire-Atlantique, France (0-m elevation)
#_______________________

# l is a list where each item corresponds to a simulatin chain GCM/RCM and contains
# a data.frame with 19 lines and 3 columns:
# - year (centered on the winter season),
# - gwl: smooth global warming levels, as shown in Fig. 2,
# - maxSWE: corresponding Annual maxima of daily SWE values (mm w.e.)
l = readRDS(file="./NUTSLoireAtlantique0m.rds")

plot(-1,-1,xlim=c(0.3,5.6),ylim=c(0,47),
     main="(a) SWE maxima in Loire-Atlantique, France (0-m elevation)",cex.main=1.5,
     xlab="Global warming level (°C)",ylab="Annual maxima of SWE (mm w.e.)", 
     cex.axis=1.3,cex.lab=1.5)

for(i in 1:length(vec.chain)){
  chain = vec.chain[i]
  df = l[[chain]]
  lines(df$gwl,df$maxSWE,lty=1,lwd=2,col=mycol.chain[[chain]])
}
grid()

legend("topright",bty="n",legend=names(mycol.chain),seg.len=4,cex=1,
       col=vec.col,lty=1,lwd=3)


#_______________________
# (b) Osttirol, Austria (2000-m elevation)
#_______________________

# l is a list where each item corresponds to a simulatin chain GCM/RCM and contains
# a data.frame with 19 lines and 3 columns:
# - year (centered on the winter season),
# - gwl: smooth global warming levels, as shown in Fig. 2,
# - maxSWE: corresponding Annual maxima of daily SWE values (mm w.e.)
l = readRDS(file="./NUTSOsttirol2000m.rds")

plot(-1,-1,xlim=c(0.3,5.6),ylim=c(0,600),
     main="(a) SWE maxima in Osttirol, Austria (2000-m elevation)",cex.main=1.5,
     xlab="Global warming level (°C)",ylab="Annual maxima of SWE (mm w.e.)", 
     cex.axis=1.3,cex.lab=1.5)

for(i in 1:length(vec.chain)){
  chain = vec.chain[i]
  df = l[[chain]]
  lines(df$gwl,df$maxSWE,lty=1,lwd=2,col=mycol.chain[[chain]])
}
grid()

dev.off()