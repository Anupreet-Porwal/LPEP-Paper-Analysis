library(ncvreg)
library(glmnet)
library(PRROC)

# Load the datasets
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/09142021_logistic_data_p100.rda?raw=TRUE"))

# Function to calculate auprc given probabilities and labels - NOT shown in paper
auprc <- function(probs, lab){
  if(all(is.na(probs))){
    return (NA)
  }else{
    probs <- probs
    fg <- probs[lab==TRUE]
    bg <- probs[lab==FALSE]
    pr <- pr.curve(scores.class0 = fg,scores.class1 = bg)

    return(pr$auc.integral)
  }
}

# Function to calculate auroc given probabilities and labels - NOT shown in paper
auroc <- function(probs, lab){
  probs <- probs
  fg <- probs[lab==TRUE]
  bg <- probs[lab==FALSE]
  roc <- roc.curve(scores.class0 = fg,scores.class1 = bg)

  return(roc$auc)
}

# Function to calculate Matthews correlation coefficient given probabilities and labels - NOT shown in paper
matthews.corr <- function(predicted, actual){
  tp <-  sum(actual==T & predicted==T)
  tn <-  sum(actual==F & predicted==F)
  fp <-  sum(actual==F & predicted==T)
  fn <-  sum(actual==T & predicted==F)

  mcc <- (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

  return(mcc)

}

# Function to calculate F1 score
F1score <- function(truegam, map.mod){
  sens <- sum(truegam ==T & map.mod==T)/sum(truegam)#sensitivity(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
  prec <- sum(truegam ==T & map.mod==T)/sum(map.mod)#posPredValue(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
  F1 <- 2*sens*prec/(sens+prec)
  return(F1)
}


set.seed(9)
nsim <- length(Beta.full)

# Matrix to store results
mcc.mat= auprc.mat =edit.dis= map.count = mse.res =avg.model.size = F1.mat <-
  matrix(0, nrow=nsim,ncol = 3)
colnames(mcc.mat)=colnames(auprc.mat)=colnames(edit.dis)= colnames(map.count) =
  colnames(mse.res) = colnames(F1.mat) = colnames(avg.model.size) <-
  c("LASSO","SCAD","MCP")



for(ind in 1:nsim){
  #### Prepare input for model runs ####
  if(ind%%100==0)cat(paste(ind," "))
  p = p.full.list[[ind]];
  betatrue = Beta.full[[ind]];
  X = X.full[[ind]][, 1:p];
  colnames(X)=paste("V",1:p,sep="")
  y = y.full[[ind]];
  sim.data <- as.data.frame(cbind(X,y))
  true.gam <- (abs(betatrue[-1])>0)
  c <- 1


  # LASSO
  lasso.fit <- cv.glmnet(X, y, family = "binomial")
  lasso.mod <- (abs(coef(lasso.fit,s="lambda.min"))>0)[-1]
  if(all(true.gam==lasso.mod)){
    map.count[ind,c] <-1
  }
  mse.res[ind,c] <- mean((betatrue-coef(lasso.fit,s="lambda.min"))^2)
  avg.model.size[ind,c] <- sum(lasso.mod)+1
  F1.mat[ind, c] <- F1score(true.gam,lasso.mod)
  edit.dis[ind,c] <- adist(toString(as.integer(lasso.mod)),toString(as.integer(true.gam)))
  auprc.mat[ind,c] <- auprc(lasso.mod,true.gam)
  mcc.mat[ind,c] <- matthews.corr(lasso.mod,true.gam)

  c=c+1
  # SCAD
  scad.fit <- cv.ncvreg(X,y,family = "binomial", penalty = "SCAD")
  scad.mod <- (abs(coef(scad.fit))>0)[-1]
  if(all(true.gam==scad.mod)){
    map.count[ind,c] <-1
  }
  mse.res[ind,c] <- mean((betatrue-coef(scad.fit))^2)
  avg.model.size[ind,c] <- sum(scad.mod)+1
  F1.mat[ind, c] <- F1score(true.gam,scad.mod)
  edit.dis[ind,c] <- adist(toString(as.integer(scad.mod)),toString(as.integer(true.gam)))
  auprc.mat[ind,c] <- auprc(scad.mod,true.gam)
  mcc.mat[ind,c] <- matthews.corr(scad.mod,true.gam)


  c=c+1
  # MCP
  mcp.fit <- cv.ncvreg(X,y,family = "binomial", penalty = "MCP")
  mcp.mod <- (abs(coef(mcp.fit))>0)[-1]
  if(all(true.gam==mcp.mod)){
    map.count[ind,c] <-1
  }
  mse.res[ind,c] <- mean((betatrue-coef(mcp.fit))^2)
  avg.model.size[ind,c] <- sum(mcp.mod)+1
  F1.mat[ind, c] <- F1score(true.gam,mcp.mod)
  edit.dis[ind,c] <- adist(toString(as.integer(mcp.mod)),toString(as.integer(true.gam)))
  auprc.mat[ind,c] <- auprc(mcp.mod,true.gam)
  mcc.mat[ind,c] <- matthews.corr(mcp.mod,true.gam)

}


# Store results as list and save them
results.all <- list("map.counts"=map.count, "mse.res"=mse.res,
                    "F1"=F1.mat, "ModelSize"=avg.model.size,
                    "edit.dis"=edit.dis,"auprc"=auprc.mat, "mcc"=mcc.mat)


save(results.all, file="./results/sim-p100/freq/alldatsets.rda")
