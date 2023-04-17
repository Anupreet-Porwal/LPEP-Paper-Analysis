rm(list=ls())

# Save working directory to address where all results are stored
setwd("./time-per-var")

# Load the datasets
load(url("https://github.com/Anupreet-Porwal/LPEP/blob/master/data/09142021_logistic_data_p100.rda?raw=TRUE"))


timepervar.summary <- function(scenario,method){
  
  s <- scenario
  r <- 1:100+(s-1)*100
  fnames <- paste(method,"/","Dataset_",r,".rda",sep = "")
  
  
  superList.time<-vector("list", length = length(r))
  
  clean_method_names <- c("LPEPE","LPEPL","LCE","LCL")
  
  for(i in 1:length(r)){
    
    load(fnames[i])
    results$time.pervar <- rbind(results$time.pervar,colSums(results$time.pervar))
    rownames(results$time.pervar)[nrow(results$time.pervar)] <- "total"
    nvar <- nrow(results$time.pervar)
    list.time <- vector("list",length=nvar)
    names(list.time) <- rownames(results$time.pervar)
    
    
    
    for (j in 1:nvar){
      list.time[[j]] <- data.frame(time=results$time.pervar[j, ], 
                                   parameters=rownames(results$time.pervar)[j],
                                   method=clean_method_names)
    }
    
    superList.time[[i]] <- do.call('rbind',list.time)
    rownames(superList.time[[i]]) <- NULL
    
  }
  
  time_table <- do.call('rbind',superList.time)
  
  
  return(time_table)
  
}

clean_method_names <- c("LPEPE","LPEPL","LCE","LCL")
time.table <- timepervar.summary(scenario = 6,method = "hyperg")

library(dplyr)

mean.time.table <-  time.table %>%
  filter(parameters !="total")%>%
  group_by(parameters,method) %>%
  summarise(mean.time=mean(time))

mean.time.table <- as.data.frame(mean.time.table)

mean.time.table$method <- factor(mean.time.table$method,levels = clean_method_names)

library(ggplot2)
ggplot(mean.time.table, aes(fill=parameters, y=mean.time, x=method)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Avg. time per iteration (in secs)")


bp <- ggplot(time.table, aes(x=method, y=time, group=method)) + 
  geom_boxplot(aes(fill=method))

bp + facet_wrap(~ parameters, ncol=2,scales="free_y")
