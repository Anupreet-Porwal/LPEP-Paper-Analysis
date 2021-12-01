# load libraries
library(caret)
library(mlbench)
library(MASS)
library(plyr)
library(haven)
library(BAS)
library(ncvreg)
library(glmnet)

# Expit function
expit <- function(x){
  return(1/(1+exp(-x)))
}

#### PEP - map model finder function ####
pep.predict.model <- function(gamsamp,estimator="HPM"){
  if(estimator=="HPM"){
    mod.sum <- count(as.data.frame(gamsamp), vars = colnames(gamsamp))
    map.mod <- mod.sum[which.max(mod.sum$freq),-ncol(mod.sum)]
  }
  return(unlist(map.mod))
}


# Predict function for test set using LPEP fit ;
# Allowed estimators are BMA, Highest probability model (HPM) and Median Probability model (MPM)
predict.lpep <- function(fit, test,estimator="BMA"){

  colnames(test)[ncol(test)] <- "Y"
  x.test <- model.matrix(Y~.,data = test)
  colnames(fit$BetaSamples) <- c("Intercept", paste("V",1:(ncol(fit$BetaSamples)-1),sep=""))

  if(estimator=="BMA"){
    predict.mat <- expit(x.test %*% t(fit$BetaSamples))
  }else if(estimator=="HPM"){
    map <- pep.predict.model(fit$GammaSamples,estimator = "HPM")
    beta.gam <- ifelse(abs(fit$BetaSamples[ ,-1])>0,1,0)
    beta.ind <- which(apply(beta.gam, 1, function(x) all.equal(x, map)) == "TRUE")
    beta.map <- fit$BetaSamples[beta.ind,]
    predict.mat <- expit(x.test %*% t(beta.map))

  }else if(estimator=="MPM"){
    mpm <- as.numeric((colMeans(fit$GammaSamples)>=0.5))
    names(mpm) <- paste("V",1:ncol(fit$GammaSamples),sep="")
    beta.gam <- ifelse(abs(fit$BetaSamples[ ,-1])>0,1,0)
    beta.ind <- which(apply(beta.gam, 1, function(x) all.equal(x, mpm)) == "TRUE")
    beta.mpm <- fit$BetaSamples[beta.ind,]
    predict.mat <- expit(x.test %*% t(beta.mpm))


  }

  return(rowMeans(predict.mat))

}


# Predict function for Fouskakis PEPs
# Estimators allowed are BMA, HPM and MPM
predict.fouskakis <- function(fit, test,estimator="BMA"){

  colnames(test)[ncol(test)] <- "Y"
  x.test <- model.matrix(Y~.,data = test)
  colnames(fit$betas) <- c("Intercept", paste("V",1:(ncol(fit$betas)-1),sep=""))

  if(estimator=="BMA"){
    predict.mat <- expit(x.test %*% t(fit$betas))
  }else if(estimator=="HPM"){
    map <- pep.predict.model(fit$gammas,estimator = "HPM")
    names(map) <- paste("V",1:ncol(fit$gammas),sep="")
    beta.gam <- ifelse(abs(fit$betas[ ,-1])>0,1,0)
    beta.ind <- which(apply(beta.gam, 1, function(x) all.equal(x, map)) == "TRUE")
    beta.map <- fit$betas[beta.ind,]
    predict.mat <- expit(x.test %*% t(beta.map))

  }else if(estimator=="MPM"){
    mpm <- as.numeric((colMeans(fit$gammas)>=0.5))
    names(mpm) <- paste("V",1:ncol(fit$gammas),sep="")
    beta.gam <- ifelse(abs(fit$betas[ ,-1])>0,1,0)
    beta.ind <- which(apply(fit$betas[ ,-1], 1, function(x) all.equal(x, mpm)) == "TRUE")
    beta.mpm <- fit$betas[beta.ind,]
    predict.mat <- expit(x.test %*% t(beta.mpm))


  }

  return(rowMeans(predict.mat))

}



# Load LPEP code
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")

# Load code for Fouskakis PEPs
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP-Paper-Analysis/blob/main/CR-DRPEP/FouskakisPEPs.R?raw=TRUE")

# Load the GUSTO dataset from Github
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/GUSTO.Rda?raw=TRUE"))



set.seed(9)


# Data reconfiguration
dataset.name <- "Gusto"

mydata <- GUSTO

x <- as.matrix(mydata[,-ncol(mydata)])

y <- as.matrix(mydata[ ,ncol(mydata)])

mydata <- cbind(x,y)
colnames(mydata)[ncol(mydata)] <- "y"
mydata <- as.data.frame(mydata)

n <- nrow(mydata)
p <- ncol(mydata)-1

#################

# It is advised to run analysis on different folds in parallel
# if you have access to cluster
# Provide the fold number and method to run as inputs to the slurm command
# The code below does the same
args = commandArgs(TRUE)

# Which method to run the analysis for ?
method <- as.numeric(args[1])

# Which fold to run analysis on ?
cv.num <- as.numeric(args[2])

# Create 10 folds
cv <- caret::createFolds(mydata$y, k=10, list=T) # Create 10 folds

fold.num <- paste("Fold",sprintf('%02d', cv.num),sep="")

# Divide data into train and test based on fold information
fold <- cv[[fold.num]]
data.train <- mydata[-fold, ] # Get the opposite of the test observations to train on
data.test <- mydata[fold, ]
x.train <- as.matrix(data.train[,-ncol(data.train)])
y.train <- data.train$y
n.train <- nrow(data.train)
y.f <- cbind(1,y.train)
x.test <- data.test[ ,-ncol(data.test)]


# List of methods grouped by type of method
method.block <- c("UIP","robust","hypergn", "crgn","drgn","crhypergn","drhypergn","freq")

# MCMC and Burn-in specification
nmc <- 131000
burn <- 10000


c=1

if(method==1){
  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("LPEP: g=n",
                    "LCL: g=n")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))

  print("----------------LPEP - g=n------------")
  # Laplace PEP -g=n
  lpep.n <- Laplace.pep(x.train,y.train,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE)
  ybma[ ,c] <- predict.lpep(lpep.n,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.lpep(lpep.n,data.test, estimator = "HPM")
  #ympm[ ,c] <- predict.lpep(lpep.n,data.test, estimator = "MPM")

  c=c+1

  # UIP
  print("----------------LCL - g=n------------")
  UIP.fit <- bas.glm( y~ ., data=data.train,method="BAS", family=binomial(link = "logit")
                      ,betaprior = g.prior(n.train))

  ybma[ ,c] <- predict(UIP.fit, x.test,estimator="BMA")$fit
  ymap[ ,c] <- predict(UIP.fit, x.test,estimator="HPM")$fit
  #ympm[ ,c] <- predict(UIP.fit, x.test,estimator="MPM")$fit

}else if(method==2){
  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("LPEP: robust",
                    "LCL: robust")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list

   # Laplace PEP - robust
  print("----------------LPEP - robust------------")
  lpep.hg <- Laplace.pep(x.train,y.train,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                         hyper="TRUE", hyper.type="robust")
  ybma[ ,c] <- predict.lpep(lpep.hg,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.lpep(lpep.hg,data.test, estimator = "HPM")
  #ympm[ ,c] <- predict.lpep(lpep.hg,data.test, estimator = "MPM")


  c=c+1

  # robust
  print("----------------LCL - robust------------")

  robust.fit <- bas.glm( y~ ., data=data.train,method="BAS", family=binomial(link = "logit")
                         ,betaprior = robust(n.train))

  ybma[ ,c] <- predict(robust.fit, x.test,estimator="BMA")$fit
  ymap[ ,c] <- predict(robust.fit, x.test,estimator="HPM")$fit
  #ympm[ ,c] <- predict(robust.fit, x.test,estimator="MPM")



}else if(method==3){
  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("LPEP: hyper-g/n",
                    "LCL: hyper-g/n")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list

  # Laplace PEP - hyper g/n
  print("----------------LPEP - hyper g/n------------")
  lpep.hgn <- Laplace.pep(x.train,y.train,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                          hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4)

  ybma[ ,c] <- predict.lpep(lpep.hgn,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.lpep(lpep.hgn,data.test, estimator = "HPM")
  #ympm[ ,c] <- predict(lpep.hgn,data.test, estimator = "MPM")


  c=c+1

  # hyper -g/n
  print("----------------LCL- hyper g/n------------")
  hypergn.fit <- bas.glm( y~ ., data=data.train,method="BAS", family=binomial(link = "logit")
                          ,betaprior = hyper.g.n(alpha=4,n.train))

  ybma[ ,c] <- predict(hypergn.fit, x.test,estimator="BMA")$fit
  ymap[ ,c] <- predict(hypergn.fit, x.test,estimator="HPM")$fit
  #ympm[ ,c] <- predict(hypergn.fit, x.test,estimator="MPM")

}else if(method==4){
  print(paste("------------",method.block[method],"----------"))


  methods.list <- c("crpep: g=n")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  #CRPEP -g=n
  cr.n <- pep.glm(y.f,x.train,iter=nmc+burn,discard=burn,family='binomial',model.prob='beta',
                  hyper = FALSE,prior.type="CRPEP")
  ybma[ ,c] <- predict.fouskakis(cr.n,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.fouskakis(cr.n,data.test, estimator = "HPM")



}else if(method==5){
  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("drpep: g=n")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  # DR PEP -g=n
  dr.n <- pep.glm(y.f, x.train,iter=nmc+burn, discard = burn, family = "binomial",
                  hyper = FALSE,prior.type="DRPEP")
  ybma[ ,c] <- predict.fouskakis(dr.n,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.fouskakis(dr.n,data.test, estimator = "HPM")


}else if(method==6){

  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("crpep: hyper-g/n")
  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  cr.hgn <- pep.glm(y.f, x.train,iter=nmc+burn, discard = burn, family = "binomial",
                    hyper = TRUE, hyper.type = "hyper-g/n",prior.type="CRPEP")
  ybma[ ,c] <- predict.fouskakis(cr.hgn,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.fouskakis(cr.hgn,data.test, estimator = "HPM")



}else if(method==7){

  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("drpep: hyper-g/n")

  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  # DR PEP - hyper g/n
  dr.hgn <- pep.glm(y.f, x.train,iter=nmc+burn, discard = burn, family = "binomial",
                    hyper = TRUE, hyper.type = "hyper-g/n",prior.type="DRPEP")
  ybma[ ,c] <- predict.fouskakis(dr.hgn,data.test, estimator = "BMA")
  ymap[ ,c] <- predict.fouskakis(dr.hgn,data.test, estimator = "HPM")


}else if(method==8){
  print(paste("------------",method.block[method],"----------"))
  methods.list <- c("lasso","scad","mcp")

  ybma <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  ymap <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  #ympm <- matrix(NA, nrow = nrow(data.test),ncol = length(methods.list))
  colnames(ybma) <- colnames(ymap) <- methods.list
  x.test.freq <- model.matrix(y~., data = data.test)
  # LASSO
  lasso.fit <- cv.glmnet(x.train, y.train, family = "binomial")
  ybma[ ,c] <- predict(lasso.fit,x.test.freq[ ,-1] , s="lambda.min",type="response")
  ymap[ ,c] <- ybma[ ,c]
  c=c+1

  # SCAD
  scad.fit <- cv.ncvreg(x.train,y.train,family = "binomial", penalty = "SCAD")
  ybma[ ,c] <- predict(scad.fit,x.test.freq[ ,-1] ,type="response")
  ymap[ ,c] <- ybma[ ,c]

  c=c+1

  #MCP
  mcp.fit <- cv.ncvreg(x.train,y.train,family = "binomial", penalty = "MCP")
  ybma[ ,c] <- predict(mcp.fit,x.test.freq[ ,-1] , type="response")
  ymap[ ,c] <- ybma[ ,c]

  c=c+1


}


# Store results as a list
res <- list("Ybma"=ybma,"Ymap"=ymap, "Ytrue"=data.test$y)

# Saved results for corresponding method(s) and fold
resultsFile = paste("./results/" ,dataset.name,"_pred/",method.block[method],"_",fold.num,".Rda",sep = "")

save(res,file = resultsFile)



