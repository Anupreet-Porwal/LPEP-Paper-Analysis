# Currently Uniform versions have been commented out; only has the Beta-binomial versions
# change the dimension of map-count (and others) and pip.mat before running uniform cases
#### Load packages ####

library(caret)      # install.packages('caret', dependencies = c("Depends", "Suggests"))
library(mlbench)
library(MASS)
library(plyr)
library(PRROC)
library(glmnet)
library(BAS)
packageVersion("BAS")
library(plyr)
library(testit)
library(devtools)
library(coda)

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
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LPEP-approx.R?raw=TRUE")

set.seed(8)

#### PEP - map model finder function #### 

pep.predict.model <- function(gamsamp,estimator="HPM"){
  if(estimator=="HPM"){
    mod.sum <- count(as.data.frame(gamsamp), vars = colnames(gamsamp))
    map.mod <- mod.sum[which.max(mod.sum$freq),-ncol(mod.sum)]
  }
  return(unlist(map.mod))
}


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

edit.dis = map.count =time.mat = mpm.count = mse.res = avg.model.size <- matrix(0, nrow=1,ncol = 4)
colnames(edit.dis) = colnames(time.mat) = colnames(map.count)= colnames(mpm.count) = 
  colnames(mse.res) = colnames(avg.model.size) <- c("LPEP-E: uip",
                                                    "LPEP-L: uip",
                                                    "LCE: uip",
                                                    "LCL: uip")
map.mat <- pip.mat <- ess.mat <- matrix(NA,nrow = p, ncol = 4)
colnames(map.mat) <- colnames(ess.mat)<- colnames(pip.mat) <- c("LPEP-E: uip",
                                                                "LPEP-L: uip",
                                                                "LCE: uip",
                                                                "LCL: uip")
time.pervar.mat <- matrix(NA, nrow=5,ncol = 4)
colnames(time.pervar.mat) <- c("LPEP-E: uip",
                               "LPEP-L: uip",
                               "LCE: uip",
                               "LCL: uip")
rownames(time.pervar.mat) <- c("gamma-delta","beta","omega","ystar","delta")


nmc=131000
burn=10000

r <- 1
c <- 1

#### Models ####
print("----------------LPEP-E BB(1,1) - uip------------")
start_time <- Sys.time()
modlpep.e <- Laplace.pep(X,y,nmc=nmc,burn=burn, 
                         model.prior = "beta-binomial", 
                         hyper="FALSE")
end_time <- Sys.time()
time.mat[1,c] <- end_time-start_time
map.lpep.e <- pep.predict.model(modlpep.e$GammaSamples,estimator = "HPM")
if(all(true.gam==map.lpep.e)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.lpep.e
pip.mat[ ,c] <- colMeans(modlpep.e$GammaSamples)
mpm.lpep.e <- (colMeans(modlpep.e$GammaSamples)>=0.5)
if(all(true.gam==mpm.lpep.e)){
  mpm.count[r,c] <-1
}
ess.mat[ ,c] <- effectiveSize(mcmc(modlpep.e$GammaSamples))
mse.res[r,c] <- mean((betatrue-colMeans(modlpep.e$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(modlpep.e$GammaSamples))
time.pervar.mat[ ,c] <- colMeans(modlpep.e$timemat)

c=c+1

print("----------------LPEP-L BB(1,1) - uip------------")
start_time <- Sys.time()
modlpep.l <- Laplace.pep.approx(X,y,nmc=nmc,burn=burn, 
                                model.prior = "beta-binomial", 
                                hyper="FALSE")

end_time <- Sys.time()
time.mat[1,c] <- end_time-start_time
map.lpep.l <- pep.predict.model(modlpep.l$GammaSamples,estimator = "HPM")
if(all(true.gam==map.lpep.l)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.lpep.l
pip.mat[ ,c] <- colMeans(modlpep.l$GammaSamples)
mpm.lpep.l <- (colMeans(modlpep.l$GammaSamples)>=0.5)
if(all(true.gam==mpm.lpep.l)){
  mpm.count[r,c] <-1
}
ess.mat[ ,c] <- effectiveSize(mcmc(modlpep.l$GammaSamples))
mse.res[r,c] <- mean((betatrue-colMeans(modlpep.l$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(modlpep.l$GammaSamples))
time.pervar.mat[ ,c] <- colMeans(modlpep.l$timemat)

c=c+1

print("----------------LCE: uip ------------")
start_time <- Sys.time()
modlce.n <- Laplace.pep(X,y,nmc=nmc,burn=burn, 
                        model.prior = "beta-binomial", 
                        hyper="FALSE",exact.mixture.g=TRUE)

end_time <- Sys.time()
time.mat[1,c] <- end_time-start_time
map.lce.n <- pep.predict.model(modlce.n$GammaSamples,estimator = "HPM")
if(all(true.gam==map.lce.n)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.lce.n
pip.mat[ ,c] <- colMeans(modlce.n$GammaSamples)
mpm.lce.n <- (colMeans(modlce.n$GammaSamples)>=0.5)
if(all(true.gam==mpm.lce.n)){
  mpm.count[r,c] <-1
}
ess.mat[ ,c] <- effectiveSize(mcmc(modlce.n$GammaSamples))
mse.res[r,c] <- mean((betatrue-colMeans(modlce.n$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(modlce.n$GammaSamples))
time.pervar.mat[ ,c] <- colMeans(modlce.n$timemat)



c=c+1

print("----------------LCL: uip ------------")
start_time <- Sys.time()
modlcl.n <- Laplace.pep.approx(X,y,nmc=nmc,burn=burn, 
                               model.prior = "beta-binomial", 
                               hyper="FALSE",exact.mixture.g=TRUE)

end_time <- Sys.time()
time.mat[1,c] <- end_time-start_time
map.lcl.n <- pep.predict.model(modlcl.n$GammaSamples,estimator = "HPM")
if(all(true.gam==map.lcl.n)){
  map.count[r,c] <-1
}
map.mat[ ,c] <- map.lcl.n
pip.mat[ ,c] <- colMeans(modlcl.n$GammaSamples)
mpm.lcl.n <- (colMeans(modlcl.n$GammaSamples)>=0.5)
if(all(true.gam==mpm.lcl.n)){
  mpm.count[r,c] <-1
}
ess.mat[ ,c] <- effectiveSize(mcmc(modlcl.n$GammaSamples))
mse.res[r,c] <- mean((betatrue-colMeans(modlcl.n$BetaSamples))^2)
avg.model.size[r,c] <- mean(rowSums(modlcl.n$GammaSamples))
time.pervar.mat[ ,c] <- colMeans(modlcl.n$timemat)


#### save results ####


results <- list("map.count"=map.count,"mpm.count"=mpm.count, "mse"=mse.res, 
                "pip.mat"=pip.mat,"model.size"=avg.model.size,
                "map.mat"=map.mat, "ess.mat"=ess.mat, "time.mat"=time.mat,
                "time.pervar"=time.pervar.mat)  


resultsFile = paste(res.folder,"Dataset_",ind,sep = "")

save(results, file = paste(resultsFile,"rda",sep = "."))

