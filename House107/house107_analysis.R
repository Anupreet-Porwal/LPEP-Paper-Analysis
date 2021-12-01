# load libraries
library(haven)
library(BAS)
library(ncvreg)
library(glmnet)
library(PRROC)
library(rms)
library(testit)




#Calculates 95% credible intervals based on posterior samples for coefficients
confint.lpep <- function(Betasamples){

  ci.bounds <- t(apply(Betasamples,2, quantile, probs=c(0.025,0.975)))
  return(ci.bounds)
}

# Load the LPEP code
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP/blob/master/R/LaplacePEP.R?raw=TRUE")

# Load code for Fouskakis PEPs
devtools::source_url("https://github.com/Anupreet-Porwal/LPEP-Paper-Analysis/blob/main/CR-DRPEP/FouskakisPEPs.R?raw=TRUE")

# Load the GUSTO dataset from Github
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/House107.Rdata?raw=TRUE"))

ind = 7
set.seed(ind)

# Data reconfiguration
mydata <- House107

x <- as.matrix(mydata[,-ncol(mydata)])

y <- as.matrix(mydata[ ,ncol(mydata)])


mydata <- cbind(x,y)
colnames(mydata)[ncol(mydata)] <- "y"
mydata <- as.data.frame(mydata)
y.f <- cbind(1,y)

n <- nrow(mydata)
p <- ncol(mydata)-1

# MCMC and Burn-in configuration
nmc <- 131000
burn <- 10000

# List of methods
# Running these methods sequentially can take a few days
# It is advised to make the code parallel and run these different
# techniques in parallel if you have access to a cluster

methods.list <- c("LPEP: g=n",
                  "LCL: g=n",
                  "CRPEP: g=n",
                  "DRPEP: g=n",
                  "LPEP: robust",
                  "LCL: robust",
                  "LPEP: hyper-g/n",
                  "LCL: hyper-g/n",
                  "CRPEP: hyper-g/n",
                  "DRPEP: hyper-g/n",
                  "LASSO","SCAD","MCP"
)

# Define matrix to store results
pip.mat <- matrix(NA,nrow=length(methods.list), ncol=p)
coef.mat <-matrix(NA,nrow=length(methods.list), ncol=p+1)
colnames(pip.mat) <- colnames(x)
colnames(coef.mat) <- c("Intercept",colnames(x))
rownames(coef.mat) <- rownames(pip.mat) <- methods.list
ci.methods <- vector("list", length(methods.list)-3)
names(ci.methods) <- methods.list[1:(length(methods.list)-3)]


c=1
# Laplace PEP -g=n
lpep.n <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper=FALSE)
pip.mat[c, ] <- colMeans(lpep.n$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.n$BetaSamples)
ci.methods[[c]] <- confint.lpep(lpep.n$BetaSamples)

c=c+1

# UIP
UIP.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                    ,betaprior = g.prior(n))
pip.mat[ c, ] <- UIP.fit$probne0[-1]
coef.mat[c, ] <- coef(UIP.fit)$postmean
ci.methods[[c]] <- confint(coef(UIP.fit))[ ,1:2]

c=c+1

#CRPEP -g=n
cr.n <- pep.glm(y.f,x,iter=nmc+burn,discard=burn,family='binomial',model.prob='beta',
                hyper = FALSE,prior.type="CRPEP")
pip.mat[c,] <- colMeans(cr.n$gammas)
coef.mat[c,  ] <- colMeans(cr.n$betas)
ci.methods[[c]] <- confint.lpep(cr.n$betas)
c=c+1

# DR PEP -g=n
dr.n <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                hyper = FALSE,prior.type="DRPEP")
pip.mat[c,] <- colMeans(dr.n$gammas)
coef.mat[c,  ] <- colMeans(dr.n$betas)
ci.methods[[c]] <- confint.lpep(dr.n$betas)
c=c+1


# Laplace PEP - robust
lpep.hg <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper="TRUE", hyper.type="robust")
pip.mat[ c, ] <- colMeans(lpep.hg$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.hg$BetaSamples)
ci.methods[[c]] <- confint.lpep(lpep.hg$BetaSamples)

c=c+1

# robust
robust.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                       ,betaprior = robust(n))
pip.mat[ c, ] <- robust.fit$probne0[-1]
coef.mat[c, ] <- coef(robust.fit)$postmean
ci.methods[[c]] <- confint(coef(robust.fit))[ ,1:2]

c=c+1

# Laplace PEP - hyper g/n
lpep.hgn <- Laplace.pep(x,y,nmc=nmc,burn=burn, model.prior = "beta-binomial", hyper="TRUE", hyper.type="hyper-g/n",hyper.param=4)
pip.mat[ c, ] <- colMeans(lpep.hgn$GammaSamples)
coef.mat[c, ] <- colMeans(lpep.hgn$BetaSamples)
ci.methods[[c]] <- confint.lpep(lpep.hgn$BetaSamples)

c=c+1

# hyper -g/n
hypergn.fit <- bas.glm( y~ ., data=mydata,method="BAS", family=binomial(link = "logit")
                        ,betaprior = hyper.g.n(alpha=4,n))
pip.mat[ c, ] <- hypergn.fit$probne0[-1]
coef.mat[c, ] <- coef(hypergn.fit)$postmean
ci.methods[[c]] <- confint(coef(hypergn.fit))[ ,1:2]

c=c+1

# CRPEP - hyper g/n
cr.hgn <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                  hyper = TRUE, hyper.type = "hyper-g/n",prior.type="CRPEP")
pip.mat[c,] <- colMeans(cr.hgn$gammas)
coef.mat[c,  ] <- colMeans(cr.hgn$betas)
ci.methods[[c]] <- confint.lpep(cr.hgn$betas)
c=c+1


# DR PEP - hyper g/n
dr.hgn <- pep.glm(y.f, x,iter=nmc+burn, discard = burn, family = "binomial",
                  hyper = TRUE, hyper.type = "hyper-g/n",prior.type="DRPEP")
pip.mat[c,] <- colMeans(dr.hgn$gammas)
coef.mat[c,  ] <- colMeans(dr.hgn$betas)
ci.methods[[c]] <- confint.lpep(dr.hgn$betas)
c=c+1


# LASSO
lasso.fit <- cv.glmnet(x, y, family = "binomial")
pip.mat[ c, ] <- (abs(coef(lasso.fit,s="lambda.min"))>0)[-1]
coef.mat[c, ] <- as.numeric(coef(lasso.fit,s="lambda.min"))
c=c+1

# SCAD
scad.fit <- cv.ncvreg(x,y,family = "binomial", penalty = "SCAD")
pip.mat[c, ] <- (abs(coef(scad.fit))>0)[-1]
coef.mat[c, ] <- coef(scad.fit)

c=c+1

#MCP
mcp.fit <- cv.ncvreg(x,y,family = "binomial", penalty = "MCP")
pip.mat[c, ] <- (abs(coef(mcp.fit))>0)[-1]
coef.mat[c, ] <- coef(mcp.fit)

c=c+1


# Store results as a list
res.house107 <-list("PIP"=pip.mat,
                "coef" = coef.mat,
                "credible.int"=ci.methods)
# Save the results
save(res.house107,file = "house107_res.Rda")

# Save the heat map plot for PIPs
pdf(file="pip_house107.pdf",height = 6,width = 8)
library(pheatmap)
pheatmap(pip.mat,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         angle_col = 45,
         color=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
