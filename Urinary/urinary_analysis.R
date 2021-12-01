# load libraries
library(BAS)
library(ncvreg)
library(glmnet)

# Calculates posterior probabilities of each model in model space
postprobs <- function(gamsamp){
    p <- ncol(gamsamp)
    model.list <-  as.matrix(expand.grid(replicate(p, 0:1, simplify = FALSE)))
    freq <- rep(0, 2^p)
    for (i in 1:nrow(model.list)){
      freq[i] <-  sum(apply(gamsamp, 1, function(x){ all(x==model.list[i, ])}))
    }
    return(freq/sum(freq))
}

# Calculates posterior probabilities for each model in model space for BAS method
bas.postprobs <- function(bas.fit){
  m <- bas.fit$n.models
  p <- bas.fit$n.vars -1
  model.probs <- bas.fit$postprobs
  models <- bas.fit$which
  model.list <-  as.matrix(expand.grid(replicate(p, 0:1, simplify = FALSE)))
  model.list.bas <- matrix(0, nrow = m, ncol = p)
  for (i in 1:m){
    if(length(models[[i]][-1])==0){next}else{
      model.list.bas[i, models[[i]][-1]] <- 1
    }
  }
  reordered.probs <- rep(0,m)
  for( i in 1:length(reordered.probs)){
    reordered.probs[i] <- model.probs[which(apply(model.list.bas, 1, function(x){ all(x==model.list[i, ])} )=="TRUE")]
  }
  return(reordered.probs)
}

#Calculates 95% credible intervals based on posterior samples for coefficients
confint.lpep <- function(Betasamples){

  ci.bounds <- t(apply(Betasamples,2, quantile, probs=c(0.025,0.975)))
  return(ci.bounds)
}

set.seed(9)
# Load LPEP code
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")
# Load the dataset
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/urinary.RData?raw=TRUE"))

# Load code for Fouskakis PEPs
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP-Paper-Analysis/blob/main/CR-DRPEP/FouskakisPEPs.R?raw=TRUE")



#Data reconfiguration
y <- urinary$y
x <- as.matrix(urinary[ ,-1])
mydata <- urinary
y.f <- cbind(1,y)

# MLE calculation - may produce warning
mle <- glm(y~.,mydata,family = binomial(link="logit"))

# List of methods
methods.list <- c("LPEP-g=n",
                  "LCL-g=n",
                  "CR PEP-g=n",
                  "DR PEP-g=n",
                  "LPEP-robust",
                  "LCL robust",
                  "LPEP-hyper-g/n",
                  "LCL hyper-g/n",
                  "CR PEP hyper-g/n",
                  "DR PEP hyper-g/n",
                  "LASSO","SCAD","MCP"
)

n <- nrow(urinary)
p <- ncol(urinary)-1

# MCMC and burn-in specifications
nmc <- 10000
burn <- 10000

# Define matrix to store results
pip.mat <- matrix(NA,nrow=length(methods.list), ncol=p)
coef.mat <-matrix(NA,nrow=length(methods.list), ncol=p+1)
colnames(pip.mat) <- paste("x",1:p, sep="")
colnames(coef.mat) <- c("Intercept",paste("x",1:p, sep=""))
rownames(coef.mat) <- rownames(pip.mat) <- methods.list
model.list <-  as.matrix(expand.grid(replicate(p, 0:1, simplify = FALSE)))
model.prob <- matrix(0,nrow = nrow(model.list),ncol = length(methods.list))
colnames(model.prob) <- rownames(coef.mat)
model.list <- cbind(model.list,model.prob)
ci.methods <- vector("list", length(methods.list)-3)
names(ci.methods) <- methods.list[1:(length(methods.list)-3)]

c=1
# Laplace PEP -g=n
lpep.n <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE)
pip.mat[c, ] <- colMeans(lpep.n$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.n$BetaSamples)
model.list[ ,c+p] <- postprobs(lpep.n$GammaSamples)
ci.methods[[c]] <- confint.lpep(lpep.n$BetaSamples)

c=c+1

# LCL UIP
UIP.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                    ,betaprior = g.prior(n))
pip.mat[ c, ] <- UIP.fit$probne0[-1]
coef.mat[c, ] <- coef(UIP.fit)$postmean
model.list[ ,c+p] <- bas.postprobs(UIP.fit)
ci.methods[[c]] <- confint(coef(UIP.fit))[ ,1:2]

c=c+1

#CRPEP -g=n
cr.n <- pep.glm(y.f,x,iter=nmc+burn,discard=burn,family='binomial',model.prob='beta',
                hyper = FALSE,prior.type="CRPEP")
pip.mat[c,] <- colMeans(cr.n$gammas)
coef.mat[c,  ] <- colMeans(cr.n$betas)
model.list[ ,c+p] <- postprobs(cr.n$gammas)
ci.methods[[c]] <- confint.lpep(cr.n$betas)
c=c+1

# DR PEP -g=n
dr.n <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                hyper = FALSE,prior.type="DRPEP")
pip.mat[c,] <- colMeans(dr.n$gammas)
coef.mat[c,  ] <- colMeans(dr.n$betas)
model.list[ ,c+p] <- postprobs(dr.n$gammas)
ci.methods[[c]] <- confint.lpep(dr.n$betas)
c=c+1



# Laplace PEP - robust
lpep.hg <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper="TRUE", hyper.type="robust")
pip.mat[ c, ] <- colMeans(lpep.hg$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.hg$BetaSamples)
model.list[ ,c+p] <- postprobs(lpep.hg$GammaSamples)
ci.methods[[c]] <- confint.lpep(lpep.hg$BetaSamples)

c=c+1

# robust
robust.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                       ,betaprior = robust(n))
pip.mat[ c, ] <- robust.fit$probne0[-1]
coef.mat[c, ] <- coef(robust.fit)$postmean
model.list[ ,c+p] <- bas.postprobs(robust.fit)
ci.methods[[c]] <- confint(coef(robust.fit))[ ,1:2]

c=c+1

# Laplace PEP - hyper g/n
lpep.hgn <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4)
pip.mat[ c, ] <- colMeans(lpep.hgn$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.hgn$BetaSamples)
model.list[ ,c+p] <- postprobs(lpep.hgn$GammaSamples)
ci.methods[[c]] <- confint.lpep(lpep.hgn$BetaSamples)

c=c+1

# hyper -g/n
hypergn.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                        ,betaprior = hyper.g.n(alpha=4,n))
pip.mat[ c, ] <- hypergn.fit$probne0[-1]
coef.mat[c, ] <- coef(hypergn.fit)$postmean
model.list[ ,c+p] <- bas.postprobs(hypergn.fit)
ci.methods[[c]] <- confint(coef(hypergn.fit))[ ,1:2]

c=c+1




# CRPEP - hyper g/n
cr.hgn <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                  hyper = TRUE, hyper.type = "hyper-g/n",prior.type="CRPEP")
pip.mat[c,] <- colMeans(cr.hgn$gammas)
coef.mat[c,  ] <- colMeans(cr.hgn$betas)
model.list[ ,c+p] <- postprobs(cr.hgn$gammas)
ci.methods[[c]] <- confint.lpep(cr.hgn$betas)
c=c+1


# DR PEP - hyper g/n
dr.hgn <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                  hyper = TRUE, hyper.type = "hyper-g/n",prior.type="DRPEP")
pip.mat[c,] <- colMeans(dr.hgn$gammas)
coef.mat[c,  ] <- colMeans(dr.hgn$betas)
model.list[ ,c+p] <- postprobs(dr.hgn$gammas)
ci.methods[[c]] <- confint.lpep(dr.hgn$betas)
c=c+1



# LASSO
lasso.fit <- cv.glmnet(x, y, family = "binomial")
pip.mat[ c, ] <- (abs(coef(lasso.fit,s="lambda.min"))>0)[-1]
coef.mat[c, ] <- as.numeric(coef(lasso.fit,s="lambda.min"))
model.selected <- which(apply(model.list[ ,1:3], 1, function(x){ all(x==pip.mat[c, ])} )=="TRUE")
model.list[model.selected,c+p ] <- 1
c=c+1

# SCAD
scad.fit <- cv.ncvreg(x,y,family = "binomial", penalty = "SCAD")
pip.mat[c, ] <- (abs(coef(scad.fit))>0)[-1]
coef.mat[c, ] <- coef(scad.fit)
model.selected <- which(apply(model.list[ ,1:3], 1, function(x){ all(x==pip.mat[c, ])} )=="TRUE")
model.list[model.selected,c+p ] <- 1

c=c+1

#MCP
mcp.fit <- cv.ncvreg(x,y,family = "binomial", penalty = "MCP")
pip.mat[c, ] <- (abs(coef(mcp.fit))>0)[-1]
coef.mat[c, ] <- coef(mcp.fit)
model.selected <- which(apply(model.list[ ,1:3], 1, function(x){ all(x==pip.mat[c, ])} )=="TRUE")
model.list[model.selected ,c+p] <- 1

c=c+1


# Store results as a list
res.urinary <- list("PIP"=pip.mat,
                    "coef" = coef.mat,
                    "model.probs"= model.list,
                    "credible.int"=ci.methods)
# Save results
save(res.urinary,file = "urinary_res.Rda")
