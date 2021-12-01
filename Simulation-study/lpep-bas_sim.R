# Currently Uniform versions have been commented out; only has the Beta-binomial versions
# change the dimension of map-count (and others) and pip.mat before running uniform cases

#### Load packages ####

library(caret)
library(mlbench)
library(MASS)
library(plyr)
library(PRROC)
library(glmnet)
library(BAS)
library(plyr)
library(testit)
library(devtools)


#### Get arguments from sbatch ####
# It is advised to run simulated datasets in parallel on cluster
# to save time

args = commandArgs(TRUE)
ind = as.numeric(args[1]) # Which simulated dataset to run analysis for
res.folder <- as.character(args[2]) # Folder where the results need to be stored
data.location <- as.character(args[3]) # Location where simulated datasets are saved

#### Load data and code ####

load(data.location)
# Load LPEP code
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")



set.seed(8)

#### PEP - map model finder function ####

pep.predict.model <- function(gamsamp,estimator="HPM"){
  if(estimator=="HPM"){
    mod.sum <- count(as.data.frame(gamsamp), vars = colnames(gamsamp))
    map.mod <- mod.sum[which.max(mod.sum$freq),-ncol(mod.sum)]
  }
  return(unlist(map.mod))
}

# Function to calculate edit distance - NOT shown in paper
edit.distance <- function(gamsamp,true.gam){
  edit.dis <- mean(apply(gamsamp, 1, function(x){adist(toString(x),toString(as.integer(true.gam)))}))
  return(edit.dis)
}

edit.distance.bas <- function(bas.fit,true.gam){
  m <- bas.fit$n.models
  p <- bas.fit$n.vars -1
  models <- bas.fit$which
  model.probs <- bas.fit$postprobs
  model.list.bas <- matrix(0, nrow = m, ncol = p)
  for (i in 1:m){
    if(length(models[[i]][-1])==0){next}else{
      model.list.bas[i, models[[i]][-1]] <- 1
    }
  }
  edit.dis <- model.probs %*% apply(model.list.bas, 1, function(x){adist(toString(x),toString(as.integer(true.gam)))})
  return(edit.dis)
}


#### Prepare input for model runs ####
p = p.full.list[[ind]];
betatrue = Beta.full[[ind]];
X = X.full[[ind]][, 1:p];
colnames(X)=paste("V",1:p,sep="")
y = y.full[[ind]];
sim.data <- as.data.frame(cbind(X,y))
true.gam <- (abs(betatrue[-1])>0)

# Define matrix to store results
edit.dis = map.count = mpm.count = mse.res =avg.model.size <- matrix(0, nrow=1,ncol = 6)
colnames(edit.dis) =  colnames(map.count)= colnames(mpm.count) = colnames(mse.res) = colnames(avg.model.size) <-
  c("LPEP-g=n","LPEP-robust",
    "LPEP-hypergn","BAS-g=n", "BAS-robust",
    "BAS-hypergn")

map.mat <- pip.mat <- matrix(NA,nrow = p, ncol = 6)
colnames(map.mat) <- colnames(pip.mat) <- c("LPEP-g=n","LPEP-robust",
                       "LPEP-hypergn","BAS-g=n", "BAS-robust",
                       "BAS-hypergn")

accept.ratios <- matrix(NA, nrow=3, ncol=3)
colnames(accept.ratios) <- c("LPEP-g=n","LPEP-robust",
                             "LPEP-hypergn")
rownames(accept.ratios) <- c("y.star.total",'y.star.local','delta')

# MCMC and burn-in specification
nmc=131000
burn=10000

r <- 1
c <- 1

#### Models ####
print("----------------LPEP- BB(1,1) - delta=n------------")
mod8 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE)
map.8 <- pep.predict.model(mod8$GammaSamples,estimator = "HPM")
if(all(true.gam==map.8)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.8
pip.mat[ ,c] <- colMeans(mod8$GammaSamples)
mpm.8 <- (colMeans(mod8$GammaSamples)>=0.5)
if(all(true.gam==mpm.8)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-colMeans(mod8$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(mod8$GammaSamples))
accept.ratios[ ,c] <- c(mod8$acc.ratio.ystar,mod8$acc.ratio.ystar.local,mod8$acc.ratio.delta)
edit.dis[r,c] <- edit.distance(mod8$GammaSamples,true.gam)

c=c+1

print("----------------LPEP- BB(1,1) - robust------------")
mod9 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                    hyper="TRUE", hyper.type="robust")
map.9 <- pep.predict.model(mod9$GammaSamples,estimator = "HPM")
if(all(true.gam==map.9)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.9
pip.mat[ ,c] <- colMeans(mod9$GammaSamples)
mpm.9 <- (colMeans(mod9$GammaSamples)>=0.5)
if(all(true.gam==mpm.9)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-colMeans(mod9$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(mod9$GammaSamples))
accept.ratios[ ,c] <- c(mod9$acc.ratio.ystar,mod9$acc.ratio.ystar.local,mod9$acc.ratio.delta)
edit.dis[r,c] <- edit.distance(mod9$GammaSamples,true.gam)


c=c+1


print("----------------LPEP- BB(1,1) - hyper-g/n------------")
mod10 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                     hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4)
map.10 <- pep.predict.model(mod10$GammaSamples,estimator = "HPM")
if(all(true.gam==map.10)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.10
pip.mat[ ,c] <- colMeans(mod10$GammaSamples)
mpm.10 <- (colMeans(mod10$GammaSamples)>=0.5)
if(all(true.gam==mpm.10)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-colMeans(mod10$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(mod10$GammaSamples))
accept.ratios[ ,c] <- c(mod10$acc.ratio.ystar,mod10$acc.ratio.ystar.local,mod10$acc.ratio.delta)
edit.dis[r,c] <- edit.distance(mod10$GammaSamples,true.gam)


c=c+1

print("----------------BAS UIP------------")

# UIP
if(p==20){
  UIP.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
                      ,betaprior = g.prior(n))

}else if(p==100){
  UIP.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
                      ,betaprior = g.prior(n),method="MCMC",
                      n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
                      laplace = FALSE)

}

UIP.hpm <- predict(UIP.fit, estimator = "HPM")
UIP.hpm.mod <- colnames(X) %in%  variable.names(UIP.hpm)
if(all(true.gam==UIP.hpm.mod)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- UIP.hpm.mod
pip.mat[ ,c] <- UIP.fit$probne0[-1]
mpm.11 <- (UIP.fit$probne0[-1]>=0.5)
if(all(true.gam==mpm.11)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-coef(UIP.fit)$postmean)^2)
avg.model.size[r,c] <- sum(UIP.fit$postprobs*UIP.fit$size)
edit.dis[r,c] <- edit.distance.bas(UIP.fit,true.gam)



c=c+1
print("----------------BAS robust------------")

# Robust
if(p==20){
  robust.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
                         ,betaprior = robust(n=n))

}else if(p==100){
  robust.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
                         ,betaprior = robust(n=n),method="MCMC",
                         n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
                         laplace = FALSE)

}


robust.hpm <- predict(robust.fit, estimator = "HPM")
robust.hpm.mod <- colnames(X) %in%  variable.names(robust.hpm)
if(all(true.gam==robust.hpm.mod)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- robust.hpm.mod
pip.mat[ ,c] <- robust.fit$probne0[-1]
mpm.12 <- (robust.fit$probne0[-1]>=0.5)
if(all(true.gam==mpm.12)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-coef(robust.fit)$postmean)^2)
avg.model.size[r,c] <- sum(robust.fit$postprobs*robust.fit$size)
edit.dis[r,c] <- edit.distance.bas(robust.fit,true.gam)



c=c+1


print("----------------BAS hyper-g/n------------")

if(p==20){
  hypergn.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
                          ,betaprior = hyper.g.n(alpha=4,n))

}else if(p==100){
  hypergn.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
                          ,betaprior = hyper.g.n(alpha=4,n),method="MCMC",
                          n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
                          laplace = FALSE)
}

hypergn.hpm <- predict(hypergn.fit, estimator = "HPM")
hypergn.hpm.mod <- colnames(X) %in%  variable.names(hypergn.hpm)
if(all(true.gam==hypergn.hpm.mod)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- hypergn.hpm.mod
pip.mat[ ,c] <- hypergn.fit$probne0[-1]
mpm.13 <- (hypergn.fit$probne0[-1]>=0.5)
if(all(true.gam==mpm.13)){
  mpm.count[r,c] <-1
}
mse.res[r,c] <- mean((betatrue-coef(hypergn.fit)$postmean)^2)
avg.model.size[r,c] <- sum(hypergn.fit$postprobs*hypergn.fit$size)
edit.dis[r,c] <- edit.distance.bas(hypergn.fit,true.gam)



c=c+1


#### save results ####


results <- list("map.count"=map.count,"mpm.count"=mpm.count, "mse"=mse.res,
                "pip.mat"=pip.mat,"model.size"=avg.model.size,
                'accept.ratio.lpep'=accept.ratios,"edit.dis"=edit.dis,
                "map.mat"=map.mat)


resultsFile = paste(res.folder,"Dataset_",ind,sep = "")

save(results, file = paste(resultsFile,"rda",sep = "."))


#### Uniform cases if we need them ####
#
# print("----------------LPEP- Uniform - delta=n------------")
# # Laplace PEP
# mod1 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform", hyper=FALSE)
#
# # colnames(mod1$GammaSamples)=paste("V",1:5,sep="")
# # # Count the unique models and their counts
# # mod.sum <- count(as.data.frame(mod1$GammaSamples), vars = paste("V",1:5,sep=""))
# # #
# # Count the unique models and their counts
# map.1 <- pep.predict.model(mod1$GammaSamples,estimator = "HPM")
# if(all(true.gam==map.1)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- colMeans(mod1$GammaSamples)
# mpm.1 <- (colMeans(mod1$GammaSamples)>=0.5)
# if(all(true.gam==mpm.1)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-colMeans(mod1$BetaSamples))^2)
# avg.model.size[r,c] <- mean(rowSums(mod1$GammaSamples))
#
#
# c=c+1
#
# print("----------------LPEP- Uniform - hyper-g------------")
# mod2 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform", hyper="TRUE", hyper.type="hyper-g",hyper.param=3)
# map.2 <- pep.predict.model(mod2$GammaSamples,estimator = "HPM")
# if(all(true.gam==map.2)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- colMeans(mod2$GammaSamples)
# mpm.2 <- (colMeans(mod2$GammaSamples)>=0.5)
# if(all(true.gam==mpm.2)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-colMeans(mod2$BetaSamples))^2)
# avg.model.size[r,c] <- mean(rowSums(mod2$GammaSamples))
#
#
# c=c+1
#
# print("---------------- LPEP- Uniform - hyper -g/n a=4------------")
# mod3 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform", hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4)
# map.3 <- pep.predict.model(mod3$GammaSamples,estimator = "HPM")
# if(all(true.gam==map.3)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- colMeans(mod3$GammaSamples)
# mpm.3 <- (colMeans(mod3$GammaSamples)>=0.5)
# if(all(true.gam==mpm.3)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-colMeans(mod3$BetaSamples))^2)
# avg.model.size[r,c] <- mean(rowSums(mod3$GammaSamples))
#
#
# c=c+1
#
#
# # UIP
# if(p==20){
#   UIP.u.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
#                         ,betaprior = g.prior(n),modelprior = uniform())
# }else if(p==100){
#   UIP.u.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
#                         ,betaprior = g.prior(n),modelprior = uniform(),method="MCMC",
#                         n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
#                         laplace = FALSE)
# }
#
# UIP.u.hpm <- predict(UIP.u.fit, estimator = "HPM")
# UIP.u.hpm.mod <- colnames(X) %in%  variable.names(UIP.u.hpm)
# if(all(true.gam==UIP.u.hpm.mod)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- UIP.u.fit$probne0[-1]
# mpm.4 <- (UIP.u.fit$probne0[-1]>=0.5)
# if(all(true.gam==mpm.4)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-coef(UIP.u.fit)$postmean)^2)
# avg.model.size[r,c] <- mean(UIP.u.fit$postprobs*UIP.u.fit$size)
#
#
# c=c+1
#
# # Hyper g
# if(p==20){
#   hyperg.u.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
#                            ,betaprior = hyper.g(alpha=3),modelprior = uniform())
# }else if(p==100){
#   hyperg.u.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
#                            ,betaprior = hyper.g(alpha=3),modelprior = uniform(),method="MCMC",
#                            n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
#                            laplace = FALSE)
#
# }
# hyperg.u.hpm <- predict(hyperg.u.fit, estimator = "HPM")
# hyperg.u.hpm.mod <- colnames(X) %in%  variable.names(hyperg.u.hpm)
# if(all(true.gam==hyperg.u.hpm.mod)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- hyperg.u.fit$probne0[-1]
# mpm.5 <- (hyperg.u.fit$probne0[-1]>=0.5)
# if(all(true.gam==mpm.5)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-coef(hyperg.u.fit)$postmean)^2)
# avg.model.size[r,c] <- mean(hyperg.u.fit$postprobs*hyperg.u.fit$size)
#
#
# c=c+1
#
# # # Hyper g/n alpha =3
# # if(p==20){
# #   hypergna3.u.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
# #                               ,betaprior = hyper.g.n(alpha=3,nrow(sim.data)),modelprior = uniform())
# #
# # }else if(p==100){
# #   hypergna3.u.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
# #                               ,betaprior = hyper.g.n(alpha=3,nrow(sim.data)),modelprior = uniform(),method="MCMC",
# #                               n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
# #                               laplace = FALSE)
# #
# # }
# #
# # hypergna3.u.hpm <- predict(hypergna3.u.fit, estimator = "HPM")
# # hypergna3.u.hpm.mod <- colnames(X) %in%  variable.names(hypergna3.u.hpm)
# # if(all(true.gam==hypergna3.u.hpm.mod)){
# #   map.count[r,c] <-1
# # }
# # c=c+1
#
# # Hyper g/n alpha =4 ; median is n
#
# if(p==20){
#   hypergn.u.fit <- bas.glm( y~ ., data=sim.data,method="BAS", family=binomial(link = "logit")
#                             ,betaprior = hyper.g.n(alpha=4,nrow(sim.data)),modelprior = uniform())
#
# }else if(p==100){
#   hypergn.u.fit <- bas.glm( y~ ., data=sim.data, family=binomial(link = "logit")
#                             ,betaprior = hyper.g.n(alpha=4,nrow(sim.data)),modelprior = uniform(),method="MCMC",
#                             n.models = nmc, MCMC.iterations = nmc, initprobs = 'eplogp',
#                             laplace = FALSE)
#
# }
#
# hypergn.u.hpm <- predict(hypergn.u.fit, estimator = "HPM")
# hypergn.u.hpm.mod <- colnames(X) %in%  variable.names(hypergn.u.hpm)
# if(all(true.gam==hypergn.u.hpm.mod)){
#   map.count[r,c] <-1
# }
# pip.mat[ ,c] <- hypergn.u.fit$probne0[-1]
# mpm.6 <- (hypergn.u.fit$probne0[-1]>=0.5)
# if(all(true.gam==mpm.6)){
#   mpm.count[r,c] <-1
# }
# mse.res[r,c] <- mean((betatrue-coef(hypergn.u.fit)$postmean)^2)
# avg.model.size[r,c] <- mean(hypergn.u.fit$postprobs*hypergn.u.fit$size)
#
# c=c+1


