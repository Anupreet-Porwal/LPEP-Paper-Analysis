rm(list=ls())

library(reshape2)
library(dplyr)

# Save working directory to address where all results are stored
setwd("./results/sim-p100/")
# setwd("./results/sim-p20/")

# Load the datasets
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/09142021_logistic_data_p100.rda?raw=TRUE"))

# Function to calculate auprc given probabilities and labels - NOT shown in paper
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

# Function to calculate auroc given probabilities and labels - NOT shown in paper
auroc <- function(probs, lab){
  probs <- probs # exclude the intercept
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


# Function to consolidate all Bayesian results
clydesim.table <- function(scenario, nmethods=13,method){
  s <- scenario
  r <- 1:100+(s-1)*100

  fnames <- paste(method,"/","Dataset_",r,".rda",sep = "")

  map.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mpm.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mse.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  edit.mat <- matrix(NA, nrow=length(r), ncol = nmethods)


  mcc.mat <- F1.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  auprc.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  modsize.mat <- matrix(NA, nrow=length(r), ncol = nmethods)

  for(i in 1:length(r)){
    load(fnames[i])
    map.mat[i, ] <- results$map.count
    mpm.mat[i, ] <- results$mpm.count
    mse.mat[i, ] <- results$mse
    edit.mat[i, ] <- results$edit.dis
    modsize.mat[i, ] <- results$model.size

    mpm.mod <- results$pip.mat>=0.5
    truegam <- abs(Beta.full[[r[i]]])>0

    for(j in 1:ncol(results$pip.mat)){
      sens <- sum(truegam[-1] ==T & mpm.mod[ ,j]==T)/sum(truegam[-1])#sensitivity(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
      prec <- sum(truegam[-1] ==T & mpm.mod[ ,j]==T)/sum(mpm.mod[ ,j])#posPredValue(as.factor(mpm.mod[ ,j]), as.factor(truegam[-1]))
      F1.mat[i,j] <- 2*sens*prec/(sens+prec)
      mcc.mat[i,j] <- matthews.corr(mpm.mod[ ,j], truegam[-1])
    }
    pip.mat <- results$pip.mat
    for(k in 1:ncol(pip.mat)){
      auprc.mat[i,k] <- auprc(pip.mat[ ,k],truegam[-1])
    }
  }

  colnames(edit.mat) <- colnames(modsize.mat) <- colnames(map.mat) <- colnames(mse.mat) <-
    colnames(mpm.mat) <- colnames(F1.mat) <- colnames(mcc.mat) <- colnames(auprc.mat) <-
    colnames(results$map.count)

  res.summary <- list("map.count"=map.mat,"mpm.count"=mpm.mat,
                      "mse"=mse.mat,"F1"=F1.mat,"auprc"=auprc.mat,"MCC"=mcc.mat,
                      "modsize"=modsize.mat, "edit.dis"=edit.mat)

  return(res.summary)
}

# Function to consolidate all frequentist results
freq.table <- function(scenario, nmethods=3,method){
  s <- scenario
  r <- 1:100+(s-1)*100

  fnames <- paste(method,"/","alldatsets",".rda",sep = "")

  map.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mpm.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mse.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  edit.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  mcc.mat <- F1.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  auprc.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  modsize.mat <- matrix(NA, nrow=length(r), ncol = nmethods)
  F1.mat <- matrix(NA, nrow=length(r), ncol = nmethods)

    load(fnames)
    map.mat <- results.all$map.count[r,]
    mpm.mat <- results.all$map.count[r,]# same as map
    mse.mat <- results.all$mse[r,]
    edit.mat <- results.all$edit.dis[r,]
    mcc.mat <- results.all$mcc[r, ]
    auprc.mat <- results.all$auprc[r, ]
    modsize.mat <- results.all$ModelSize[r, ]
    F1.mat <- results.all$F1[r, ]

    colnames(F1.mat) <- colnames(modsize.mat) <- colnames(map.mat) <- colnames(mse.mat) <-
      colnames(mpm.mat) <-colnames(F1.mat) <- colnames(mcc.mat) <- colnames(auprc.mat) <-
      colnames(results.all$map.count)

  res.summary <- list("map.count"=map.mat,"mpm.count"=mpm.mat,
                      "mse"=mse.mat,"F1"=F1.mat,"auprc"=auprc.mat,"MCC"=mcc.mat,
                      "modsize"=modsize.mat,"edit.dis"=edit.mat)

  return(res.summary)
}


method_folders <- c("lpep-bas","exact-g","freq")
nmethods <- c(6,3,3)

# Matrix to store consolidated results
auprc.mat <-  modsize.mat <- edit.mat <- mcc.mat <- F1.mat <- map.count <- mse.mat <- mpm.count  <-
  matrix(NA, nrow = sum(nmethods), ncol = 8)
colnames(auprc.mat) <-colnames(modsize.mat) <-colnames(edit.mat) <- colnames(mcc.mat) <- colnames(F1.mat) <-
  colnames(map.count) <- colnames(mse.mat) <- colnames(mpm.count) <-
  paste("scenario",1:8, sep = "")


datalist=list()
p=list()
for (j in 1:ncol(map.count)){
    # For each scenario; consolidate results over 100 datasets
    res.summary.lpep <- clydesim.table(scenario = j,nmethods = nmethods[1],method=method_folders[1])
    res.summary.exact <- clydesim.table(scenario = j,nmethods = nmethods[2],method=method_folders[2])
    res.summary.freq <- freq.table(scenario = j,nmethods = nmethods[3],method=method_folders[3])

    # Calculate summary
    map.count[ ,j] <- c(colSums(res.summary.lpep$map.count),colSums(res.summary.exact$map.count),colSums(res.summary.freq$map.count))
    mpm.count[ ,j] <- c(colSums(res.summary.lpep$mpm.count), colSums(res.summary.exact$mpm.count),colSums(res.summary.freq$mpm.count))
    mse.mat[ ,j] <- c(colMeans(res.summary.lpep$mse),colMeans(res.summary.exact$mse),colMeans(res.summary.freq$mse))
    F1.mat[, j] <- c(colMeans(res.summary.lpep$F1),colMeans(res.summary.exact$F1),colMeans(res.summary.freq$F1))
    mcc.mat[ ,j] <- c(colMeans(res.summary.lpep$MCC),colMeans(res.summary.exact$MCC),colMeans(res.summary.freq$MCC))
    edit.mat[ ,j] <- c(colMeans(res.summary.lpep$edit.dis),colMeans(res.summary.exact$edit.dis),colMeans(res.summary.freq$edit.dis))
    modsize.mat[ ,j] <- c(colMeans(res.summary.lpep$modsize),colMeans(res.summary.exact$modsize),colMeans(res.summary.freq$modsize))
    auprc.mat[ ,j] <- c(colMeans(res.summary.lpep$auprc),colMeans(res.summary.exact$auprc),colMeans(res.summary.freq$auprc))

}



rownames(edit.mat) <- rownames(mcc.mat) <- rownames(F1.mat) <- rownames(map.count) <-
  rownames(mse.mat) <- rownames(mpm.count) <- rownames(auprc.mat) <- rownames(modsize.mat) <-
  c(colnames(res.summary.lpep$map.count),colnames(res.summary.exact$map.count),colnames(res.summary.freq$map.count))

# Reorder methods to have order similar to paper
order <- c(1,4,7,2,5,8,3,6,9,10,11,12)

# Tables corresponding to paper tables
auprc.mat <- auprc.mat[order, ]
F1.mat <- F1.mat[order, ]
mcc.mat <- mcc.mat[order, ]
map.count <- map.count[order, ]
mse.mat <- mse.mat[order, ]
mpm.count <- mpm.count[order, ]
edit.mat <- edit.mat[order, ]
modsize.mat <- modsize.mat[order, ]


###### F1 plot
datalist=list()
for (j in 1:8){

  if(j %in% c(1,2)){ptrue=0}
  if(j %in% c(3,4)){ptrue=5}
  if(j %in% c(5,6)){ptrue=10}
  if(j %in% c(7,8)){ptrue=20}
  res.summary.lpep <- clydesim.table(scenario = j,nmethods = nmethods[1],method=method_folders[1])
  res.summary.exact <- clydesim.table(scenario = j,nmethods = nmethods[2],method=method_folders[2])
  res.summary.freq <- freq.table(scenario = j,nmethods = nmethods[3],method=method_folders[3])

  if(j>=3){
    F1.agg.uip <- cbind(res.summary.lpep$F1[ ,1],res.summary.exact$F1[ ,1],res.summary.lpep$F1[ ,4])
    F1.agg.robust <- cbind(res.summary.lpep$F1[ ,2],res.summary.exact$F1[ ,2],res.summary.lpep$F1[ ,5])
    F1.agg.hypergn <- cbind(res.summary.lpep$F1[ ,3],res.summary.exact$F1[ ,3],res.summary.lpep$F1[ ,6])
    colnames(F1.agg.uip) <- colnames(F1.agg.robust) <- colnames(F1.agg.hypergn) <- c('LPEP','LCE', 'LCL')
    F1.agg.freq <- res.summary.freq$F1

    m1 <- cbind.data.frame(label=rep('UIP',100),F1.agg.uip)
    m2 <- cbind.data.frame(label=rep('Robust',100),F1.agg.robust)
    m3 <- cbind.data.frame(label=rep('Hyper-g/n',100),F1.agg.hypergn)
    m4 <- cbind.data.frame(label=rep('Frequentist',100),F1.agg.freq)

    df <- rbind.data.frame(m1,m2,m3)

    df.m <- melt(df, id.var = "label")
    df.mfreq <- melt(m4,id.vars = "label")

    df.final <- rbind(df.m,df.mfreq)
    df.final$label <- factor(df.final$label,levels =
                               c("UIP","Robust","Hyper-g/n","Frequentist"))
    datalist[[j-2]] <- cbind.data.frame(df.final,
                                        corr=rep(ifelse(j %% 2==1,"r=0","r=0.75" ),100),
                                        ptrue=rep(ptrue,100))

  }

}

big_data = do.call(rbind, datalist)
ptrue.labs <- c("p-true= 5","p-true= 10","p-true= 20")
names(ptrue.labs) <- c('5','10','20')

# pdf(file="F1plot_p20.pdf",height = 12,width = 8)
# ggplot(data = big_data, aes(x=label, y=value,fill=variable)) +
#   geom_boxplot(outlier.shape = 21,outlier.size = 1) +
#   stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red",
#                position = position_dodge2(width = 0.75,
#                                           preserve = "single"))+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size = 5))+
#   xlab("method")+ylab("F1 score")+
#   ylim(0,1) + facet_grid(ptrue~corr,labeller = labeller(ptrue=ptrue.labs))
# dev.off()

pdf(file="F1plot.pdf",height = 12,width = 8)
ggplot(data = big_data, aes(x=label, y=value,fill=variable)) +
  geom_boxplot(outlier.shape = 21,outlier.size = 1) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red",
               position = position_dodge2(width = 0.75,
                                          preserve = "single"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  xlab("method")+ylab("F1 score")+
  ylim(0,1) + facet_grid(ptrue~corr,labeller = labeller(ptrue=ptrue.labs))
dev.off()




###### Model size plot
datalist=list()
for (j in 1:8){

  if(j %in% c(1,2)){ptrue=0}
  if(j %in% c(3,4)){ptrue=5}
  if(j %in% c(5,6)){ptrue=10}
  if(j %in% c(7,8)){ptrue=20}
  res.summary.lpep <- clydesim.table(scenario = j,nmethods = nmethods[1],method=method_folders[1])
  res.summary.exact <- clydesim.table(scenario = j,nmethods = nmethods[2],method=method_folders[2])
  res.summary.freq <- freq.table(scenario = j,nmethods = nmethods[3],method=method_folders[3])

  modsize.agg.uip <- cbind(res.summary.lpep$modsize[ ,1],res.summary.exact$modsize[ ,1],res.summary.lpep$modsize[ ,4])
  modsize.agg.robust <- cbind(res.summary.lpep$modsize[ ,2],res.summary.exact$modsize[ ,2],res.summary.lpep$modsize[ ,5])
  modsize.agg.hypergn <- cbind(res.summary.lpep$modsize[ ,3],res.summary.exact$modsize[ ,3],res.summary.lpep$modsize[ ,6])
  colnames(modsize.agg.uip) <- colnames(modsize.agg.robust) <- colnames(modsize.agg.hypergn) <- c('LPEP','LCE', 'LCL')
  modsize.agg.freq <- res.summary.freq$modsize

  m1 <- cbind.data.frame(label=rep('UIP',100),modsize.agg.uip)
  m2 <- cbind.data.frame(label=rep('Robust',100),modsize.agg.robust)
  m3 <- cbind.data.frame(label=rep('Hyper-g/n',100),modsize.agg.hypergn)
  m4 <- cbind.data.frame(label=rep('Frequentist',100),modsize.agg.freq)

  df <- rbind.data.frame(m1,m2,m3)

  df.m <- melt(df, id.var = "label")
  df.mfreq <- melt(m4,id.vars = "label")

  df.final <- rbind(df.m,df.mfreq)
  df.final$label <- factor(df.final$label,levels =
                             c("UIP","Robust","Hyper-g/n","Frequentist"))
  datalist[[j]] <- cbind.data.frame(df.final,
                                    corr=rep(ifelse(j %% 2==1,"r=0","r=0.75" ),100),
                                    ptrue=rep(ptrue,100))
}

big_data = do.call(rbind, datalist)
ptrue.labs <- c("p-true= 0","p-true= 5","p-true= 10","p-true= 20")
names(ptrue.labs) <- c('0','5','10','20')

pdf(file="modsizeplot.pdf",height = 12,width = 8)
ggplot(data = big_data, aes(x=label, y=value,fill=variable)) +
  geom_boxplot(outlier.shape = 21,outlier.size = 1) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red",
               position = position_dodge2(width = 0.75,
                                          preserve = "single"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  xlab("method")+ylab("Model Size")+
  ylim(0,100) +
  facet_grid(ptrue~corr,labeller = labeller(ptrue=ptrue.labs))+
  geom_hline(aes(yintercept = as.numeric(ptrue)),linetype = "dashed",color="blue")
dev.off()

# # For simulation where p =20 (Supplementary materials);
# #results can be requested from the authors

# pdf(file="modsizeplot_p20.pdf",height = 12,width = 8)
# ggplot(data = big_data, aes(x=label, y=value,fill=variable)) +
#   geom_boxplot(outlier.shape = 21,outlier.size = 1) +
#   stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red",
#                position = position_dodge2(width = 0.75,
#                                           preserve = "single"))+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size = 5))+
#   xlab("method")+ylab("Model Size")+
#   ylim(0,20) +
#   facet_grid(ptrue~corr,labeller = labeller(ptrue=ptrue.labs))+
#   geom_hline(aes(yintercept = as.numeric(ptrue)),linetype = "dashed",color="blue")
# dev.off()


##### MSE plot

datalist=list()
for (j in 1:8){

  if(j %in% c(1,2)){ptrue=0}
  if(j %in% c(3,4)){ptrue=5}
  if(j %in% c(5,6)){ptrue=10}
  if(j %in% c(7,8)){ptrue=20}
  res.summary.lpep <- clydesim.table(scenario = j,nmethods = nmethods[1],method=method_folders[1])
  res.summary.exact <- clydesim.table(scenario = j,nmethods = nmethods[2],method=method_folders[2])
  res.summary.freq <- freq.table(scenario = j,nmethods = nmethods[3],method=method_folders[3])

  mse.agg.uip <- cbind(res.summary.lpep$mse[ ,1],res.summary.exact$mse[ ,1],res.summary.lpep$mse[ ,4])*1000
  mse.agg.robust <- cbind(res.summary.lpep$mse[ ,2],res.summary.exact$mse[ ,2],res.summary.lpep$mse[ ,5])*1000
  mse.agg.hypergn <- cbind(res.summary.lpep$mse[ ,3],res.summary.exact$mse[ ,3],res.summary.lpep$mse[ ,6])*1000
  colnames(mse.agg.uip) <- colnames(mse.agg.robust) <- colnames(mse.agg.hypergn) <- c('LPEP','LCE', 'LCL')
  mse.agg.freq <- res.summary.freq$mse*1000

  m1 <- cbind.data.frame(label=rep('UIP',100),mse.agg.uip)
  m2 <- cbind.data.frame(label=rep('Robust',100),mse.agg.robust)
  m3 <- cbind.data.frame(label=rep('Hyper-g/n',100),mse.agg.hypergn)
  m4 <- cbind.data.frame(label=rep('Frequentist',100),mse.agg.freq)

  df <- rbind.data.frame(m1,m2,m3)

  df.m <- melt(df, id.var = "label")
  df.mfreq <- melt(m4,id.vars = "label")

  df.final <- rbind(df.m,df.mfreq)
  df.final$label <- factor(df.final$label,levels =
                             c("UIP","Robust","Hyper-g/n","Frequentist"))
  datalist[[j]] <- cbind.data.frame(df.final,
                                    corr=rep(ifelse(j %% 2==1,"r=0","r=0.75" ),100),
                                    ptrue=rep(ptrue,100))


}

big_data = do.call(rbind, datalist)
ptrue.labs <- c("p-true=0","p-true= 5","p-true= 10","p-true= 20")
names(ptrue.labs) <- c('0','5','10','20')

filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]

  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

big_data <- big_data %>% mutate(big_data, method=paste(label,variable,sep="+"))
pdf(file="mseplot.pdf",height = 12,width = 8)
big_data %>%
  group_by(corr,ptrue,method)%>%
  mutate(value2 = filter_lims(value)) %>%
  ggplot(aes(x=label, y=value2,fill=variable)) +
  geom_boxplot(na.rm = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red",
               position = position_dodge2(width = 0.75,
                                          preserve = "single"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  xlab("method")+ylab("1000*AMSE")+
  facet_grid(ptrue~corr,labeller = labeller(ptrue=ptrue.labs),scales = "free")
dev.off()

