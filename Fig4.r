# ---
# title: "Figure 4: Illustration of return levels for the NUTS-3 
# Loire-Atlantique, France. Each column corresponds to a nonstationary model,
# and each row to a simulation chain. Return levels are shown for different 
# annual exceedance probabilities as a function of the global warming level. 
# Red boxes indicate the preferred model. For each simulation chain, annual SWE
# maxima are shown with black points."
# author: "Guillaume Evin"
# date: "02/12/2024"
# ---

# This script produces the Figure 4 in the article entitled "Estimating changes 
# in extreme snow load conditionally to global warming levels" by Evin. G.,
# Le Roux, E., Kamir, E. and Morin, S., submitted to Cold Regions Science and 
# Technology.

# lib_fitSWE.r contains the functions that fit the nonstationary models, calculate
# the return levele estimates, etc.
source("./lib_fitSWE.r")

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
n.chain = length(vec.chain)

# threshold for the SWE (mm)
threshold = 1

# l is a list where each item corresponds to a simulatin chain GCM/RCM and contains
# a data.frame with 19 lines and 3 columns:
# - year (centered on the winter season),
# - gwl: smooth global warming levels, as shown in Fig. 2,
# - maxSWE: corresponding Annual maxima of daily SWE values (mm w.e.)
l = readRDS(file="./NUTSLoireAtlantique0m.rds")

# we compute the return levels for warming levels ranging from 0?C to 6?C
pred.gt = seq(from=0,to=6,by=0.1)

# vector of return levels and corresponding colors
pred.RL = c(0.5,0.75,0.9,0.98)
col.RL = c("#1b9e77","#d95f02","#7570b3","#e7298a")
n.RL = length(pred.RL)

# for this NUTS3 area, three models are compared: pEXP, pGA and pIG
n.model = 3

pdf("../Figures/Fig4_exFitWithZeros.pdf",width = 8,height = 8)
par(mfrow=c(n.chain,n.model),mar=c(0,0,0,0),oma=c(4,4,2,10))
for(i.chain in 1:n.chain){
  # GCM / RCM
  chain = vec.chain[i.chain]
  
  # reformat the data.frame and consider values smaller than the threshold as
  # zero
  df = data.frame(gt=l[[chain]]$gwl,maxswe=l[[chain]]$maxSWE)
  df.zz = get.df.nona(df,threshold)
  
  # fit.SWE fit all the candidate model and return the outputs of the
  # gamlss and gamlssx packages
  lmodel = fit.SWE(df.zz)
  
  # i.best.model select the model with the minimum AIC
  i.best = i.best.model(lmodel)
  
  for(i.model in 1:n.model){
    selmod = lmodel[[i.model]]
    
    # first plot the annual SWE maxima (black points)
    plot(df$gt,df$maxswe,pch=20,col="black",
         xlab="",ylab="",ylim=c(0,30),xaxt="n",yaxt="n",frame.plot = F)
    if(i.model==1){
      legend("topleft",legend = chain,cex=0.7,bty="n")
    }
    grid()
    
    # add return level curves
    for(i.RL in 1:n.RL){
      # get.predicted.rl returns the return level estimate corresponding to 
      # a probability "proba" and a vector of global warming levels
      RL = get.predicted.rl(selmod,pred.gt=pred.gt,proba=pred.RL[i.RL],df.zz)
      lines(pred.gt,RL,col=col.RL[i.RL],lwd=2)
    }
  
    # add a red box if it is the selected model
    if(i.model==i.best){
      box(col = 'red',lwd=3)
    }
    
    # add axis
    if(i.chain==n.chain){
      axis(side=1,at = seq(from=0,to=6,by=1), tck = 0.05, padj=-1)
    }
    if(i.model==1){
      axis(side=2,at = seq(from=10,to=30,by=10), tck = 0.05, padj=1)
    }
  }
  
}
mtext(text = c("pEXP","pGA","pIG"),
      side = 3, at=seq(from=0.15,to=0.85,length.out=n.model), outer = T,line = 0.5)
mtext(text = "Global warming levels (°C)",side = 1, outer = T,line = 2.5,cex=1.5)
mtext(text = "Return levels of SWE (mm w.e.)",side = 2, outer = T,line = 2,cex=1.5)


par(fig = c(0.75, 1, 0.45, 0.8), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottomright",bty="n",legend=c("0.50","0.75","0.90","0.98"),
       title="Probability",
       seg.len=4,cex=1.2,col=col.RL,lty=1,lwd=3)

dev.off()