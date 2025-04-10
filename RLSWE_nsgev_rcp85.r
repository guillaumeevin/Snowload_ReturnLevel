# Guillaume Evin
# 09/04/2025
# guillaume.evin@inrae.fr
#
# fit non-stationary models to maxima of SWE
#
# - if there are zero-values, we fit mixed gamma distributions (p + (1-p)f_gamma)
# - if there are zero-values, we fit gamma distributions, gumbels or gev distributions
#
# - AIC criteria will be used to select the best model
# - different fit for each model couples GCM/RCM
# - global temperature as predictors
# - mix the different scenarios
# - we can consider linear or smoothing splines for the evolution of the par.
# - 27/08/2024: add p-value of a statistic test

library(parallel) #detectCores
library(foreach) #parallelization
library(doParallel) #parallelization
library(doSNOW)
library(boot)

#===============================================================================
# load data
#===============================================================================

#Number of cores for parallelization
nbcore=detectCores()-2 #Number of cores for parallelization
cl <- makeSOCKcluster(nbcore)
registerDoSNOW(cl)

# period for the models (global temperatures available until 2099)
vy = 1951:2099

maxSWE = readRDS("../MTMSI-max-swe/maxSWE.rds")
dfRuns = readRDS("dfRuns.rds")
source("lib_fitSWE.r")
n.scen = dim(maxSWE)[1]
n.nutsalt = dim(maxSWE)[2]
sl.years = 1951:2100 # checked
zz.sl = sl.years%in%vy


# load global temperature: vecYearsGT,HadCRUT5.spline,GCM.spline,globalTempRef,scen_GCM_RCP,scenAvail,gtas
load("../GlobalTAS/GlobalTempSmooth.RData")
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


#===============================================================================
# fit models on the nuts and get 50-year return levels
#===============================================================================
pb <- txtProgressBar(min=1, max=n.nutsalt, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# get couple GCM/RCM
gcm.rcm = paste0(dfRuns$gcm,"/",dfRuns$rcm)
couple.gcm.rcm = unique(gcm.rcm)
n.couple = length(couple.gcm.rcm)

# write rsd file for each nut to avoid referencing the big matrix in foreach
for(inut in 1:n.nutsalt){
   saveRDS(object = maxSWE[,inut,zz.sl],file = paste0("../MMENUTS/maxSWE",inut,".rds"))
}

# Creating a function to pass into boot() function
bootFunc <- function(mydata, i, labelBestModel){
  df <- mydata[i,]
  tryCatch(
    expr = {
      bootmodel = fit.1nsmodel(data = df,label = labelBestModel)
      get.predicted.rl(bootmodel,pred.gt=c(1,2,3),proba=0.98,data = df)
    },
    error = function(cond){
      NA
    },
    warning = function(cond){
      NA
    }
  )   
}

fitmodel <- function(data){
  tryCatch(
    expr = {
      fit.SWE(data = data)
    },
    error = function(cond){
      NA
    },
    warning = function(cond){
      NA
    }
  )    
}

# loop over the nut
# inut=4017 # 2000m/Osttirol/Austria
# inut=4935 # 0m/Loire-Atlantique
threshold = 1 # https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/2016GL071789
l=foreach(inut=1:n.nutsalt, .packages = c('evd','gamlss','gamlssx','boot'), .options.snow=opts) %dopar% {
  maxSWEnut = readRDS(file = paste0("../MMENUTS/maxSWE",inut,".rds")) # 20 scenarios * 149 years
  rl = rl2nd = p0 = matrix(nrow = n.couple, ncol=3)
  vec.ibest = vec.i2nd = vec.pval = vec.diffAIC = vec.shape.GEV = vector(length=n.couple)
  listmodel = list()
  for(i.couple in 1:n.couple){
    # get ony the scenario RCP85 for that GCM/RCM
    is = which(gcm.rcm==couple.gcm.rcm[i.couple]&dfRuns$rcp=="RCP85")
    gt = as.vector(matGT[,is])
    maxSWE.couple = maxSWEnut[is,]
    
    # vectorize
    maxswe = as.vector(t(maxSWE.couple))
    
    # we fit model if we do not have more than 90% of zero SWE max
    # we consider a threshold of 1 mm
    if(mean(maxSWEnut<threshold,na.rm=T)<0.9){
      df = data.frame(gt=gt,maxswe=maxswe)
      
      # filter na and set low values to 0
      df.zz = get.df.nona(df,threshold)
      
      
      # apply all ns-gev models and return a list containing models' outputs
      lmodel = fitmodel(df.zz)
      
      if(is.na(lmodel)[1]){
        returnNA = TRUE
      }else{
        returnNA = FALSE
        
        # select best model according to the AIC criteria (index in the list)
        nmodel = length(lmodel)
        vec.AIC = vector(length=nmodel)
        for(i.m in 1:nmodel){
          vec.AIC[i.m] = lmodel[[i.m]]$aic
        }
        
        i.best = which.min(vec.AIC)
        i.2nd = order(vec.AIC)[2]
        vec.diffAIC[i.couple] = vec.AIC[i.2nd]-vec.AIC[i.best]
        
        
        # return return levels, p0, label of the best model
        selmodel = lmodel[[i.best]]
        p0[i.couple,] = get.predicted.p0(model.out = selmodel,pred.gt=c(1,2,3),data = df.zz)
        rl[i.couple,] = get.predicted.rl(model.out = selmodel,pred.gt=c(1,2,3),proba=0.98,data=df.zz)
        vec.ibest[i.couple] = selmodel$label
        
        vec.i2nd[i.couple] = lmodel[[i.2nd]]$label
        rl2nd[i.couple,] = get.predicted.rl(model.out = lmodel[[i.2nd]],
                                            pred.gt=c(1,2,3),proba=0.98,data=df.zz)
        
        # test normality of the residuals
        res = residuals(selmodel)
        shapiro.out = shapiro.test(res)
        vec.pval[i.couple] = shapiro.out$p.value
        
        # add shape par GEV
        if(selmodel$label=="GEV"){
          vec.shape.GEV[i.couple] = selmodel$nu.coefficients
        }
      }
    }else{
      returnNA = TRUE
    }
    
    if(returnNA){
      vec.shape.GEV[i.couple] = NA
      p0[i.couple,] = NA
      rl[i.couple,] = NA
      vec.ibest[i.couple] = 0
      vec.i2nd[i.couple] = 0
      rl2nd[i.couple,] = NA
      vec.diffAIC[i.couple] = NA
      vec.pval[i.couple] = NA
    }
  }
  l = list(ibest=vec.ibest,pval=vec.pval,rl=rl,p0=p0,vec.shape.GEV=vec.shape.GEV,
           vec.diffAIC=vec.diffAIC,rl2nd=rl2nd,vec.i2nd=vec.i2nd)
  saveRDS(object = l,file = paste0("../ReturnLevelOutput/l",inut,".rds"))
  
  # bootstrap of return level estimats with the best model
  nboot = 100
  bootout = array(dim=c(n.couple,nboot,3))
  for(i.couple in 1:n.couple){
    if(vec.ibest[i.couple]!=0){
      is = which(gcm.rcm==couple.gcm.rcm[i.couple]&dfRuns$rcp=="RCP85")
      gt = as.vector(matGT[,is])
      maxSWE.couple = maxSWEnut[is,]
      maxswe = as.vector(t(maxSWE.couple))
      if(mean(maxSWEnut<threshold,na.rm=T)<0.9){
        df = data.frame(gt=gt,maxswe=maxswe)
        df.zz = get.df.nona(df,threshold)
        b <- boot(df.zz, bootFunc, R = nboot, labelBestModel = vec.ibest[i.couple])
        bootout[i.couple,,]=b$t
      }else{
        bootout[i.couple,,] = NA
      }
    }
  }
  b = list(bootout=bootout)
  saveRDS(object = b,file = paste0("../ReturnLevelOutput/b",inut,".rds"))
  
  NULL
}
