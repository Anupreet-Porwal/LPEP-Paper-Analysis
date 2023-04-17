# Load libraries
library(AUC);
library(PRROC)

# Function to calculate AUROC given probabilities/score x and label l
auroc <- function(x,l){ auc(roc(x,l))}


# Function to calculate Calibration slope given probabilities/score x and label l
cs <- function(x,l){glm(l ~ log(x / (1 -x)), family = binomial(link = logit))$coef[2]}

# Function to calculate Log score given probabilities/score x and label l
ls <- function(x,l){-mean(l * log(x) + (1 - l) * log(1 - x))}

# Function to calculate Brier score given probabilities/score x and label l
bs <- function(x,l){mean((l - x)^2)}


# Set working directory to where the Gusto prediction results are stored
setwd("./results/Gusto_pred/")

all.files <- list.files()

# List of methods by type
method.block <- c("UIP","robust","hypergn",  "crgn","drgn",
                  "crhypergn",
                  "drhypergn","freq","LPEPL-UIP","LPEPL-robust","LPEPL-hypergn")

# List of methods
methods.list <- c("LPEP: g=n",
                  "LCL: g=n",
                  "LPEP: robust",
                  "LCL: robust",
                  "LPEP: hyper-g/n",
                  "LCL: hyper-g/n",
                  "CRPEP: g=n",
                  "DRPEP: g=n",
                  "CRPEP: hyper-g/n",
                  "DRPEP: hyper-g/n",
                  "LASSO","SCAD","MCP",
                  "LPEPL: g=n",
                  "LPEPL: robust",
                  "LPEPL: hyper-g/n"
)

folds <- 10

fold.num <- paste("Fold",sprintf('%02d', 1:folds),sep="")

# Matrix for storing calibration slope, log score, Brier score and AUC for MAP/BMA models
 cs.map <- cs.bma <- ls.map <- ls.bma <- bs.map <- bs.bma <-  auc.map <- auc.bma <- 
   auprc.map <- auprc.bma <- 
  matrix(NA, nrow = length(methods.list),ncol = length(fold.num))

 rownames(cs.map) <- rownames(cs.bma) <- rownames(ls.map) <- rownames(ls.bma) <-
  rownames(bs.map) <- rownames(bs.bma) <-rownames(auc.map) <- rownames(auc.bma) <- methods.list
colnames(auprc.bma) <- colnames(auprc.map) <-colnames(cs.map) <- colnames(cs.bma) <- colnames(ls.map) <- colnames(ls.bma) <-
  colnames(bs.map) <- colnames(bs.bma) <-colnames(auc.map) <- colnames(auc.bma) <- fold.num


# Consolidate results from all methods
count=1
for (i in 1: length(method.block)){

  method.files <- all.files[grep(paste("^",method.block[i],sep=""), all.files)]


  for(j in 1:folds){

    if(!file.exists(paste(method.block[i],"_",fold.num[j],".Rda",sep=""))){
      next
    }else{
      load(paste(method.block[i],"_",fold.num[j],".Rda",sep=""))


      r <- count:(count+ncol(res$Ybma)-1)
      # AUC
      auc.bma[r, j] <- apply(res$Ybma, 2,auroc,l=factor(res$Ytrue))
      auc.map[r, j] <- apply(res$Ymap, 2,auroc,l=factor(res$Ytrue))

      # CS
      cs.bma[r, j] <- apply(res$Ybma,2,cs, l=res$Ytrue)
      cs.map[r, j] <- apply(res$Ymap,2,cs, l=res$Ytrue)

      # LS
      ls.bma[r,j] <- apply(res$Ybma,2,ls,l=res$Ytrue)
      ls.map[r,j] <- apply(res$Ymap,2,ls,l=res$Ytrue)

      # Brier score
      bs.bma[r,j] <- apply(res$Ybma,2,bs,l=res$Ytrue)
      bs.map[r,j] <- apply(res$Ymap,2,bs,l=res$Ytrue)

    }
  }

  count=count+ncol(res$Ybma)
}

# Create summary tables for BMA based results
summ.pred <- cbind(round(rowMeans(auc.bma,na.rm=TRUE),digits=4),
                   round(rowMeans(cs.bma,na.rm=TRUE),digits=4),
                   round(rowMeans(ls.bma,na.rm=TRUE),digits=4),
                   round(rowMeans(bs.bma,na.rm=TRUE),digits=4))
colnames(summ.pred) <- c("AUC","CS","LS","Brier")

# Create summary tables for MAP based results - NOT shown in paper
summ.pred.map <- cbind(round(rowMeans(auc.map,na.rm=TRUE),digits=4),
                   round(rowMeans(cs.map,na.rm=TRUE),digits=4),
                   round(rowMeans(ls.map,na.rm=TRUE),digits=4),
                   round(rowMeans(bs.map,na.rm=TRUE),digits=4))
colnames(summ.pred.map) <- c("AUC","CS","LS","Brier")



