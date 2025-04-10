# ---
# title: "Figure 9: Illustration of the uncertainty associated with 50-year 
# return level estimates as a function of the global warming level, for the 
# illustrative NUTS-3 areas: Loire-Atlantique, France (0 m) on the left column
# and Osttirol, Austria (2000 m) on the right column. The first nine lines 
# correspond to the different simulation chains. The different thin curves show 
# the 100 bootstrap estimates and the thick curve shows the estimate obtained
# using all SWE maxima (i.e. not bootstrapped). The last line shows the overall
# estimates. The interval in grey corresponds to the quantiles 0.05 and 0.95 
# obtained from the $100 \times 9$ estimates of 50-year return levels represented
# by all the thin curves of the first nine lines. The thick black curve is the
# mean return level estimate, averaged over the different simulation chains. For
# each simulation chain, annual SWE maxima are shown with black points.
# author: "Guillaume Evin"
# date: "02/12/2024"
# ---

# This script produces the Figure 9 in the article entitled "Estimating changes 
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


# boot is the package used for the bootstrap
library(boot)

# Creating a function to pass into boot() function
nboot = 100
bootFunc <- function(mydata, i, labelBestModel){
  df <- mydata[i,]
  tryCatch(
    expr = {
      bootmodel = fit.1nsmodel(data = df,label = labelBestModel)
      get.predicted.rl(bootmodel,pred.gt=pred.gt,proba=0.98,data = df)
    },
    error = function(cond){
      NA
    },
    warning = function(cond){
      NA
    }
  )   
}

#_______________________________________________________________________________
#  get return level estimates with bootstrap estimates
#_______________________________________________________________________________

# two NUTS3 
vec.nut = c("LoireAtlantique0m","Osttirol2000m")

# initiliaze lists that will contain the results
rl.boot = list()
rl.hat = list()

# vector of global warming levels
pred.gt = seq(from=0,to=6,by=0.1)
ngt = length(pred.gt)

for(i in 1:length(vec.nut)){
  # l is a list where each item corresponds to a simulatin chain GCM/RCM and contains
  # a data.frame with 19 lines and 3 columns:
  # - year (centered on the winter season),
  # - gwl: smooth global warming levels, as shown in Fig. 2,
  # - maxSWE: corresponding Annual maxima of daily SWE values (mm w.e.)
  l = readRDS(file=paste0("./NUTS",vec.nut[i],".rds"))

  # prepare outputs
  rl.boot[[i]] = array(dim=c(n.chain,nboot,ngt))
  rl.hat[[i]] = matrix(nrow = n.chain,ncol = ngt)
  
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
    
    # select the best model
    selmod = lmodel[[i.best]]
    
    # get 50-year return level estimates
    rl.hat[[i]][i.chain,] = get.predicted.rl(selmod,pred.gt=pred.gt,proba=0.98,df.zz)
    
    # bootstrap of return level estimats with the best model
    b <- boot(df.zz, bootFunc, R = nboot, labelBestModel = selmod$label)
    rl.boot[[i]][i.chain,,] = b$t
  }
}


#_______________________________________________________________________________
# figure: one column per NUTS3 area, one line per simulation chain
# the last line gives the overall estimate and the overall uncertainty [5%, 95%]
#_______________________________________________________________________________
pdf("../Figures/Fig9_exUnc.pdf",width = 8,height = 10)
par(mfcol=c(n.chain+1,2),mar=c(0,3,0.5,0),oma=c(12,3,2,0.5))
for(i in 1:2){
  # y limits
  if(i==1){
    myylim = c(0,40)
  }else{
    myylim = c(0,800)
  }
  
  # l is a list where each item corresponds to a simulatin chain GCM/RCM and contains
  # a data.frame with 19 lines and 3 columns:
  # - year (centered on the winter season),
  # - gwl: smooth global warming levels, as shown in Fig. 2,
  # - maxSWE: corresponding Annual maxima of daily SWE values (mm w.e.)
  l = readRDS(file=paste0("./NUTS",vec.nut[i],".rds"))
  
  for(i.chain in 1:n.chain){
    # GCM / RCM
    chain = vec.chain[i.chain]
    gt = l[[chain]]$gwl
    maxSWE = l[[chain]]$maxSWE
    
    # reformat the data.frame and consider values smaller than the threshold as
    # zero
    df = data.frame(gt=l[[chain]]$gwl,maxswe=l[[chain]]$maxSWE)
    
    # colors
    colchain = mycol.chain[[vec.chain[i.chain]]]
    colshade = adjustcolor(colchain,alpha.f=0.1)
    
    # retrieve return levels
    rl = rl.hat[[i]][i.chain,]
    rlb = rl.boot[[i]][i.chain,,]
    
    # plot
    plot(-1,-1,xlab="",ylab="",xlim=c(0,5),xaxt="n",ylim=myylim,
         frame.plot = F)
    grid()
    
    # get ony the scenario RCP85 for that GCM/RCM
    points(gt,maxSWE,pch=20,cex=1)
    
    # add the curves corresponding to all th bootstrap estimates
    for(i.boot in 1:nboot){
      lines(pred.gt,rlb[i.boot,],col=colshade,lwd=1)
    }
    
    # best estimate
    lines(pred.gt,rl,col=colchain,lwd=2)
    if(i==1){
      legend("topleft",bty="n",legend=vec.chain[i.chain])
    }
  }
  
  # all estimate
  plot(-1,-1,xlab="",ylab="",xlim=c(0,5),ylim=myylim,frame.plot = F)
  grid()
  
  # interval bootstrap
  binf = apply(rl.boot[[i]],3,quantile,probs = 0.05,na.rm=T)
  bsup = apply(rl.boot[[i]],3,quantile,probs = 0.95,na.rm=T)
  polygon(x = c(pred.gt,rev(pred.gt)),y = c(binf,rev(bsup)),col="grey",border = NA)
  
  # mean estimate
  lines(pred.gt,apply(rl.hat[[i]],2,mean),col="black",lwd=3)
  
  if(i==1){
    legend("topleft",bty="n",legend="Overall estimate")
  }
}

mtext(text = c("Loire-Atlantique, France (0 m)","Osttirol, Austria (2000 m)"),
      side = 3, at=c(0.28,0.75), outer = T,line = 0.5)
mtext(text = "Global warming levels (?C)",side = 1, outer = T,line = 3,cex=1.5)
mtext(text = "Return levels of SWE (mm w.e.)",side = 2, outer = T,line = 1,cex=1.5)

# legend
par(fig = c(0.2, 0.8, 0, 0.2), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",bty="n",
       legend=c("Annual SWE maxima","Best estimate","Bootstrap estimates","90% uncertainty interval"),
       title="",seg.len=4,cex=1.5,col=c("black","black","gray90","gray"),
       lwd=c(NA,3,1,10),lty=c(NA,1,1,1),pch=c(20,NA,NA,NA))
dev.off()