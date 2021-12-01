library(BAS)
library(ncvreg)
library(glmnet)

postprobs <- function(gamsamp){
    p <- ncol(gamsamp)
    model.list <-  as.matrix(expand.grid(replicate(p, 0:1, simplify = FALSE)))
    freq <- rep(0, 2^p)
    for (i in 1:nrow(model.list)){
      freq[i] <-  sum(apply(gamsamp, 1, function(x){ all(x==model.list[i, ])}))
    }
    return(freq/sum(freq))
}

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

confint.lpep <- function(Betasamples){

  ci.bounds <- t(apply(Betasamples,2, quantile, probs=c(0.025,0.975)))
  return(ci.bounds)
}

set.seed(9)

devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/urinary.RData?raw=TRUE"))

y <- urinary$y
x <- as.matrix(urinary[ ,-1])
mydata <- urinary

mle <- glm(y~.,mydata,family = binomial(link="logit"))

methods.list <- c("LPEP-g=n",
                  "UIP-g=n",
                  "LPEP-robust",
                  "robust",
                  #"LPEP-hyperg",
                  #"hyperg",
                  "LPEP-hyper-g/n",
                  "hyper-g/n",
                  "CR PEP hyper-g",
                  "DR PEP hyper-g",
                  "CR PEP hyper-g/n",
                  "DR PEP hyper-g/n",
                  "CR PEP-g=n",
                  "DR PEP-g=n",
                  "SH-g=n",
                  "SH-robust",
                  #"SH-hyperg",
                  "SH-hyper-g/n",
                  "LASSO","SCAD","MCP"
)#,
#"Ex-g=n",
#"Ex-hyper-g",
#"Ex-hyper-g/n")


n <- nrow(urinary)
p <- ncol(urinary)-1


nmc <- 10000
burn <- 10000

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

# UIP
UIP.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                    ,betaprior = g.prior(n))
pip.mat[ c, ] <- UIP.fit$probne0[-1]
coef.mat[c, ] <- coef(UIP.fit)$postmean
model.list[ ,c+p] <- bas.postprobs(UIP.fit)
ci.methods[[c]] <- confint(coef(UIP.fit))[ ,1:2]

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


# Fouskakis PEPs
source("C:/Users/Anupreet Porwal/Dropbox/Research/PEP-GLM/code/main/FouskakisPEPs.R")

y.f <- cbind(1,y)

# CRPEP - hyperg
cr.hg <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                 hyper = TRUE, hyper.type = "hyper-g",prior.type="CRPEP")
pip.mat[c,] <- colMeans(cr.hg$gammas)
coef.mat[c,  ] <- colMeans(cr.hg$betas)
model.list[ ,c+p] <- postprobs(cr.hg$gammas)
ci.methods[[c]] <- confint.lpep(cr.hg$betas)
c=c+1

# DR PEP - hyper g
dr.hg <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                 hyper = TRUE, hyper.type = "hyper-g",prior.type="DRPEP")
pip.mat[c,] <- colMeans(dr.hg$gammas)
coef.mat[c,  ] <- colMeans(dr.hg$betas)
model.list[ ,c+p] <- postprobs(dr.hg$gammas)
ci.methods[[c]] <- confint.lpep(dr.hg$betas)
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



# Sebanes Held -g=n
sebheld.n <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE,seb.held=TRUE)
pip.mat[c, ] <- colMeans(sebheld.n$GammaSamples)
coef.mat[c, ] <- colMeans(sebheld.n$BetaSamples)
model.list[ ,c+p] <- postprobs(sebheld.n$GammaSamples)
ci.methods[[c]] <- confint.lpep(sebheld.n$BetaSamples)
c=c+1

# Sebanes Held- robust
sebheld.robust <- Laplace.pep(x,y,nmc=nmc,burn=burn,
                          model.prior = "beta-binomial", hyper="TRUE",
                          hyper.type="robust",seb.held = TRUE)
pip.mat[ c, ] <- colMeans(sebheld.robust$GammaSamples)
coef.mat[c, ] <- colMeans(sebheld.robust$BetaSamples)
model.list[ ,c+p] <- postprobs(sebheld.robust$GammaSamples)
ci.methods[[c]] <- confint.lpep(sebheld.robust$BetaSamples)
c=c+1

# Sebanes Held- hyper g/n
sebheld.hgn <- Laplace.pep(x,y,nmc=nmc,burn=burn,
                           model.prior = "beta-binomial", hyper="TRUE",
                           hyper.type="hyper-g/n",hyper.param=4,seb.held = TRUE)
pip.mat[ c, ] <- colMeans(sebheld.hgn$GammaSamples)
coef.mat[c, ] <- colMeans(sebheld.hgn$BetaSamples)
model.list[ ,c+p] <- postprobs(sebheld.hgn$GammaSamples)
ci.methods[[c]] <- confint.lpep(sebheld.hgn$BetaSamples)
c=c+1

# # Exact BAS -g=n
# ex.n <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE,exact.mixture.g =TRUE)
# pip.mat[c, ] <- colMeans(ex.n$GammaSamples)
# coef.mat[c, ] <- colMeans(ex.n$BetaSamples)
# model.list[ ,c+p] <- postprobs(ex.n$GammaSamples)
# c=c+1
#
# # Exact BAS - hyper g
# ex.hg <- Laplace.pep(x,y,nmc=nmc,burn=burn,
#                      model.prior = "beta-binomial", hyper="TRUE",
#                      hyper.type="hyper-g",hyper.param=3,exact.mixture.g =TRUE)
# pip.mat[ c, ] <- colMeans(ex.hg$GammaSamples)
# coef.mat[c, ] <- colMeans(ex.hg$BetaSamples)
# model.list[ ,c+p] <- postprobs(ex.hg$GammaSamples)
# c=c+1
#
# # Exact BAS - hyper g/n
# ex.hgn <- Laplace.pep(x,y,nmc=nmc,burn=burn,
#                       model.prior = "beta-binomial", hyper="TRUE",
#                       hyper.type="hyper-g/n",hyper.param=4,exact.mixture.g =TRUE)
# pip.mat[ c, ] <- colMeans(ex.hgn$GammaSamples)
# coef.mat[c, ] <- colMeans(ex.hgn$BetaSamples)
# model.list[ ,c+p] <- postprobs(ex.hgn$GammaSamples)

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


order2 <- c(1,2,11,12,13,3,4,7,8,14,5,6,9,10,15,16,17,18)
pip.mat.ord <- pip.mat[order2, ]
coef.mat.ord <- coef.mat[order2,]
model.list.ord <- model.list[ ,c(1:p,order2+p)]
ci.methods <- ci.methods[order2[1:(length(order2)-3)]]

res.urinary <- list("PIP"=pip.mat.ord,"coef" = coef.mat.ord,"model.probs"= model.list.ord, "credible.int"=ci.methods)

save(res.urinary,file = "C:/Users/Anupreet Porwal/Dropbox/Research/PEP-GLM/code/results/urinary_res.Rda")
