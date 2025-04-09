library(gamlss)
library(gamlssx)

################################################################################
# get.predicted.p0
#
# 
# INPUTS: 
# - model.out: object gamlss or gamlssx
# - pred.gt: future gt
# - data: data.frame with gt and maxswe
get.predicted.p0 = function(model.out,pred.gt,data){
  if(model.out$label%in%c("pIG","pGA","pEXP")){
    p0 <- predict(model.out,what="nu", data=data, 
                  newdata = data.frame(gt = pred.gt),type="response")
  }else{
    p0 = rep(0,length(pred.gt))
  }
  
  return(p0)
}

################################################################################
# get.predicted.rl
#
# 
# INPUTS: 
# - model.out: object gamlss
# - pred.gt: future gt
# - data: data.frame with gt and maxswe
# - proba: proba for return level
get.predicted.rl = function(model.out,pred.gt,data,proba=0.98){
  # get predicted values for mu
  pred.mu <- predict(model.out,what="mu", data=data, 
                     newdata = data.frame(gt = pred.gt),type="response")
  # get predicted values for nu
  if("nu"%in%model.out$parameters){
    pred.nu <- predict(model.out,what="nu", data=data, 
                       newdata = data.frame(gt = pred.gt),type="response")
  }else{
    pred.nu = NULL
  }
  
  # get predicted values for sigma
  if("sigma"%in%model.out$parameters){
    if(!is.null(model.out$sigma.formula)){
      pred.sigma <- predict(model.out,what="sigma", data=data, 
                            newdata = data.frame(gt = pred.gt),type="response")
    }else{
      pred.sigma = rep(model.out$sigma.fv[1],length(pred.gt))
    }
  }else{
    pred.sigma = 1
  }
  
  # get return levels using the quantile functions of each distribution
  dist.family = model.out$family[1]
  q.out = vector(length=length(pred.gt))
  for(igt in 1:length(pred.gt)){
    q.out[igt] = switch(dist.family,
                   GA = qGA(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt]),
                   IG = qIG(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt]),
                   GG = qGG(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt],nu = pred.nu[igt]),
                   LOGNO = qLOGNO(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt]),
                   PARETO2 = qPARETO2(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt]),
                   ZAGA = qZAGA(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt],nu = pred.nu[igt]),
                   ZAIG = qZAIG(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt],nu = pred.nu[igt]),
                   RG = qRG(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt]),
                   GEV = qGEV(p = proba,mu = pred.mu[igt],sigma = pred.sigma[igt],nu = pred.nu[igt]))
  }
  return(q.out)
}



################################################################################
# i.best.model
#
# 
# INPUTS: 
# - list.models
i.best.model = function(list.models){
  n.models = length(list.models)
  vec.AIC = vector(length=n.models)
  for(i.m in 1:n.models){
    vec.AIC[i.m] = list.models[[i.m]]$aic
  }
  
  return(which.min(vec.AIC))
}




################################################################################
# get.df.nona
#
# filter NA values and set values below a threshold to 0
#
# INPUTS: 
# - df: data.frame with gt and maxswe data
# - threshold: threshold (in mm) below which the SWE maximum is considered as zero
get.df.nona = function(df,threshold=1){
  # filter NAs
  zz = !is.na(df$maxswe)
  df.zz = df[zz,]
  
  # set small values to zero
  df.zz$maxswe[df.zz$maxswe<threshold] = 0
  
  return(df.zz)
}


################################################################################
# fit.1nsmodel
#
# fit one non-stationary distribution
#
# INPUTS: 
# - data: data.frame with gt and maxswe data
# - label: string indicating the type of distribution to be fitted
fit.1nsmodel = function(data,label){
  if(label=="pGA"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,nu.fo=~gt,family=ZAGA,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="pEXP"){
    m = gamlss(maxswe~gt,mu.fo=~gt,nu.fo=~gt,family=ZAGA,
               data=data,sigma.fix = T,sigma.start = 1,method=mixed(50,500), trace=FALSE)
  }else if(label=="pIG"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,nu.fo=~gt,family=ZAIG,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="GA"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,family=GA,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="GG"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,family=GG,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="LOGNO"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,family=LOGNO,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="Pareto2 linear"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,family=PARETO2,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="GUM"){
    m = gamlss(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,family=RG,
               data=data,method=mixed(50,500), trace=FALSE)
  }else if(label=="GEV"){
    m = fitGEV(maxswe~gt,mu.fo=~gt,sigma.fo=~gt,
               data=data,method=mixed(50,500), trace = FALSE)
  }
  
  m$label = label
  return(m)
}

################################################################################
# fit.SWE
#
# fit non-stationary models to SWE maxima as a function of global warming
#
# CASE 1: for chains with zero values, mixed p-gamma or mixed p-exp; linear
# CASE 2: for chains with no zero values, gamma,gumbel,gev (xi cst); linear
# 
# INPUTS: 
# - data is a data.frame with two columns: maxswe and gt (smooth global temperature)
fit.SWE = function(data){
  
  # are there zero values
  hasZero = any(data$maxswe==0)
  
  list.models = list()
  if(hasZero){
    vec.lab = c("pEXP","pGA","pIG")
    for(i.lab in 1:length(vec.lab)){
      list.models[[i.lab]] = fit.1nsmodel(data = data,label = vec.lab[i.lab])
    }
  }else{
    vec.lab = c("GA","GG","GUM","GEV")
    for(i.lab in 1:length(vec.lab)){
      list.models[[i.lab]] = fit.1nsmodel(data = data,label = vec.lab[i.lab])
    }
  }
  
  return(list.models)
}