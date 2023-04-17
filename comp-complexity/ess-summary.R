rm(list=ls())

# Save working directory to address where all results are stored
rho <- 0.9
setwd("./sim-ess")

# Load the datasets
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/09142021_logistic_data_p100.rda?raw=TRUE"))

auprc <- function(probs, lab){
  if(all(is.na(probs))){
    return (NA)
  }else{
    probs <- probs # exclude the intercept
    fg <- probs[lab==TRUE]
    bg <- probs[lab==FALSE]
    pr <- pr.curve(scores.class0 = fg,scores.class1 = bg)
    
    return(pr$auc.integral)
  }
}

auroc <- function(probs, lab){
  probs <- probs # exclude the intercept
  fg <- probs[lab==TRUE]
  bg <- probs[lab==FALSE]
  roc <- roc.curve(scores.class0 = fg,scores.class1 = bg)
  
  return(roc$auc)
}


matthews.corr <- function(predicted, actual){
  tp <-  sum(actual==T & predicted==T)
  tn <-  sum(actual==F & predicted==F)
  fp <-  sum(actual==F & predicted==T)
  fn <-  sum(actual==T & predicted==F)
  
  mcc <- (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  return(mcc)
  
}


method_folders <- c("uip","hyperg")
nmethods <- c(4,4)


ess.summary <- function(scenario,method,nmethods=4){
  
  s <- scenario
  r <- 1:100+(s-1)*100
  fnames <- paste(method,"/","Dataset_",r,".rda",sep = "")
  
  map.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mpm.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mse.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  modsize.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  F1.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  
  median.ess.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  median.esr.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  time.mat <- matrix(NA, nrow = length(r),ncol=nmethods)
  median.ess.nzero <- median.esr.nzero <- matrix(NA, nrow=length(r), ncol = nmethods)
  median.ess.zero <- median.esr.zero <-  matrix(NA, nrow=length(r), ncol = nmethods)
  
  max.ess.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  max.esr.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  
  for(i in 1:length(r)){
    
    load(fnames[i])
    map.mat[i, ] <- results$map.count
    mpm.mat[i, ] <- results$mpm.count
    mse.mat[i, ] <- results$mse
    modsize.mat[i, ] <- results$model.size
    median.ess.mat[i, ] <- apply(results$ess.mat, 2, median)
    max.ess.mat[i, ] <- apply(results$ess.mat, 2, max)
    esr <- apply(results$ess.mat,1,function(x){x*3600/results$time.mat})
    time.mat[i, ] <- results$time.mat
    median.esr.mat[i, ] <- t(apply(esr, 1, median))
    max.esr.mat[i, ] <- t(apply(esr, 1, max))
    
    median.ess.nzero[i, ] <- apply(results$ess.mat[c(1:5,11:15),], 2, median)
    median.ess.zero[i, ] <- apply(results$ess.mat[-c(1:5,11:15),], 2, median)
    median.esr.nzero[i, ] <- apply(esr[,c(1:5,11:15)], 1, median)
    median.esr.zero[i, ] <- apply(esr[,-c(1:5,11:15)], 1, median)
    
    mpm.mod <- results$pip.mat>=0.5
    truegam <- abs(Beta.full[[r[i]]])>0
    
    for(j in 1:ncol(results$pip.mat)){
      sens <- sum(truegam[-1] ==T & mpm.mod[ ,j]==T)/sum(truegam[-1])#sensitivity(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
      prec <- sum(truegam[-1] ==T & mpm.mod[ ,j]==T)/sum(mpm.mod[ ,j])#posPredValue(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
      F1.mat[i,j] <- 2*sens*prec/(sens+prec) 
      #mcc.mat[i,j] <- matthews.corr(mpm.mod[ ,j], truegam[-1])
    }
  }
  
  colnames(modsize.mat) <- colnames(map.mat) <- colnames(mse.mat) <- 
    colnames(mpm.mat) <- colnames(median.ess.mat) <- colnames(max.ess.mat) <- 
    colnames(median.esr.mat) <- colnames(max.esr.mat) <- colnames(F1.mat) <- 
    colnames(median.esr.nzero) <- colnames(median.ess.nzero) <- 
    colnames(median.ess.zero) <-  colnames(median.esr.zero) <- colnames(time.mat) <- 
    c("LPEP-E","LPEP-L","LCE","LCL")
  
  res.summary <- list("map.count"=map.mat,"mpm.count"=mpm.mat,
                      "mse"=mse.mat,
                      "F1"=F1.mat,
                      "time"=time.mat,
                      "modsize"=modsize.mat, 
                      "median.ess"=median.ess.mat,
                      "max.ess"=max.ess.mat,
                      "median.esr"=median.esr.mat,
                      "max.esr"=max.esr.mat,
                      "ess.zero"=median.ess.zero,
                      "ess.nzero"=median.ess.nzero,
                      "esr.zero"=median.esr.zero,
                      "esr.nzero"=median.esr.nzero)
  
  return(res.summary)
  
}


res.summary.uip <- ess.summary(scenario=6, method = method_folders[1])
res.summary.hyperg <- ess.summary(6, method = method_folders[2])


results.hyperg <- rbind(apply(res.summary.hyperg$median.ess,2,median), 
      apply(res.summary.hyperg$median.esr,2,median),
      apply(res.summary.hyperg$time,2,median),
      apply(res.summary.hyperg$mse,2,median)*1000,
      apply(res.summary.hyperg$F1,2,median))
rownames(results.hyperg) <- c("ESS","ESR","time","MSE","F1")
#colnames(results.hyperg) <- c("LPEP-E","LPEP-L","LCE","LCL")

results.uip <- rbind(apply(res.summary.uip$median.ess,2,median), 
      apply(res.summary.uip$median.esr,2,median),
      apply(res.summary.uip$time,2,median),
      apply(res.summary.uip$mse,2,median)*1000,
      apply(res.summary.uip$F1,2,median))
rownames(results.uip) <- c("ESS","ESR","time","MSE","F1")
#colnames(results.uip) <- c("LPEP-E","LPEP-L","LCE","LCL")


# code to generate box plots
library(reshape)
ess.uip <- melt(as.data.frame(res.summary.uip$median.ess))
esr.uip <- melt(as.data.frame(res.summary.uip$median.esr))
time.uip <- melt(as.data.frame(res.summary.uip$time/3600))
f1.uip <- melt(as.data.frame(res.summary.uip$F1))
mse.uip <- melt(as.data.frame(res.summary.uip$mse))


ess.hyperg <- melt(as.data.frame(res.summary.hyperg$median.ess))
esr.hyperg <- melt(as.data.frame(res.summary.hyperg$median.esr))
time.hyperg <- melt(as.data.frame(res.summary.hyperg$time/3600))
f1.hyperg <- melt(as.data.frame(res.summary.hyperg$F1))
mse.hyperg <- melt(as.data.frame(res.summary.hyperg$mse))

par(mfrow=c(2,3))
boxplot(data=ess.uip, value~variable,ylab = "ESS",xlab = "method",outline=FALSE)
boxplot(data=time.uip, value~variable,ylab = "Time",xlab = "method",outline=FALSE)
boxplot(data=esr.uip, value~variable,ylab = "ESR",xlab = "method",outline=FALSE)
boxplot(data=f1.uip, value~variable,ylab = "F1 score",xlab = "method",outline=FALSE)
boxplot(data=mse.uip, value~variable,ylab = "MSE",xlab = "method",outline=FALSE)


par(mfrow=c(2,3))
boxplot(data=ess.hyperg, value~variable,ylab = "ESS",xlab = "method",outline=FALSE)
boxplot(data=time.hyperg, value~variable,ylab = "Time",xlab = "method",outline=FALSE)
boxplot(data=esr.hyperg, value~variable,ylab = "ESR",xlab = "method",outline=FALSE)
boxplot(data=f1.hyperg, value~variable,ylab = "F1 score",xlab = "method",outline=FALSE)
boxplot(data=mse.hyperg, value~variable,ylab = "MSE",xlab = "method",outline=FALSE)





