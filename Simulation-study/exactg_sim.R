# Load the libraries
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

devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")

set.seed(8)


#samp <- 100

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

#### PEP - map model finder function ####

pep.predict.model <- function(gamsamp,estimator="HPM"){
  if(estimator=="HPM"){
    mod.sum <- count(as.data.frame(gamsamp), vars = colnames(gamsamp))
    map.mod <- mod.sum[which.max(mod.sum$freq),-ncol(mod.sum)]
  }
  return(unlist(map.mod))
}

# Define matrix to store results
edit.dis=map.count = mpm.count = mse.res = avg.model.size <- matrix(0, nrow=1,ncol = 3)
colnames(edit.dis)= colnames(map.count)= colnames(mpm.count) =
  colnames(mse.res)= colnames(avg.model.size) <- c("Exact-g=n", "Exact-robust","Exact-hypergn")

#### Prepare input for model runs ####
p = p.full.list[[ind]];
betatrue = Beta.full[[ind]];
X = X.full[[ind]][, 1:p];
colnames(X)=paste("V",1:p,sep="")
y = y.full[[ind]];
sim.data <- as.data.frame(cbind(X,y))
true.gam <- (abs(betatrue[-1])>0)

map.mat <- pip.mat <- matrix(NA,nrow = p, ncol = 3)
colnames(map.mat) <- colnames(pip.mat) <-c("Exact-g=n", "Exact-robust","Exact-hypergn")

# MCMC and Burn-in specification
nmc <- 131000
burn <- 10000

r <- 1
c <- 1

##### Methods
print("----------------Exact- BB(1,1) - delta=n------------")
mod8 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                    hyper=FALSE,exact.mixture.g=TRUE)
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
edit.dis[r,c] <- edit.distance(mod8$GammaSamples,true.gam)


c=c+1

print("----------------Exact- BB(1,1) - robust -----------")
mod9 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                    hyper="TRUE", hyper.type="robust",exact.mixture.g=TRUE)
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
edit.dis[r,c] <- edit.distance(mod9$GammaSamples,true.gam)


c=c+1


print("----------------Exact- BB(1,1) - hyper-g/n------------")
mod10 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "beta-binomial",
                     hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4,exact.mixture.g=TRUE)
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
edit.dis[r,c] <- edit.distance(mod10$GammaSamples,true.gam)


c=c+1

#### save results ####

results <- list("map.count"=map.count,"mpm.count"=mpm.count, "mse"=mse.res,
                "pip.mat"=pip.mat,"model.size"=avg.model.size,"edit.dis"=edit.dis,
                "map.mat"=map.mat)

resultsFile = paste(res.folder,"Dataset_",ind,sep = "")

filename=paste(resultsFile, "rda",sep=".")


save(results, file = paste(resultsFile,"rda",sep = "."))


#### Code for Uniform model prior
# print("---------------- Exact Uniform - delta=n------------")
# # Laplace PEP
# mod1 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform", hyper=FALSE,
#                     exact.mixture.g=TRUE)
#
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
# c=c+1
#
# print("----------------Exact- Uniform - hyper-g------------")
# mod2 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform", hyper="TRUE",
#                     hyper.type="hyper-g",hyper.param=3,exact.mixture.g=TRUE)
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
# print("---------------- Exact- Uniform - hyper -g/n a=4------------")
# mod3 <- Laplace.pep(X,y,nmc=nmc,burn=burn, model.prior = "Uniform",
#                     hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4,exact.mixture.g=TRUE)
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
