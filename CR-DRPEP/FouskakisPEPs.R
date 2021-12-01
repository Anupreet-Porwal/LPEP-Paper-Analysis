library(mvtnorm)

expit <- function(x){
  return(1/(1+exp(-x)))
}

b.accept.prob <- function(beta.new,pi.new,beta.old,pi.old,success,trials,
                          ini.y,trials.star,
                          Beta.prop,Sigma.prop,inv.delta){

  gam1 <- length(beta.new)
  log.num <- sum(dbinom(success,trials,pi.new,log=TRUE))+
    inv.delta*sum(dbinom(ini.y,trials.star,pi.new,log=TRUE))-
    dmvnorm(beta.new,Beta.prop,Sigma.prop,log=TRUE)
  log.den <- sum(dbinom(success,trials,pi.old,log=TRUE))+
    inv.delta*sum(dbinom(ini.y,trials.star,pi.old,log=TRUE))-
    dmvnorm(beta.old,Beta.prop,Sigma.prop,log=TRUE)

  a.p = min(exp(log.num-log.den),1)

  return (a.p)
}

b0.accept.prob <- function(beta0.new,beta0.old,ini.y,trials.star,Beta0.prop,
                           Sigma0.prop,psi){
  n.star <- length(trials.star)
  pi0.new <- expit(beta0.new)
  pi0.old <- expit(beta0.old)
  pi0.new <- rep(pi0.new,length.out=n.star)
  pi0.old <- rep(pi0.old,length.out=n.star)
  log.num<-1/psi*sum(dbinom(ini.y,trials.star,pi0.new,log=TRUE))-
    dnorm(beta0.new,Beta0.prop,Sigma0.prop,log=TRUE)
  log.den<-1/psi*sum(dbinom(ini.y,trials.star,pi0.old,log=TRUE))-
    dnorm(beta0.old,Beta0.prop,Sigma0.prop,log=TRUE)
  a.p = min(exp(log.num-log.den),1)

  return (a.p)

}

ystar.accept.prob <- function(ystar.new,ystar.old, trials_star,
                              p.gamma, p0, piML.new,
                              piML.old, Det.new, Det.old, p.prop,psi,inv.delta){
  n.star <- length(trials_star)
  p.gamma<-rep(p.gamma,length.out=n.star)
  p0<-rep(p0,length.out=n.star)

  log.num <- inv.delta*sum(dbinom(ystar.new,trials_star,p.gamma,log=TRUE))+
    1/psi*sum(dbinom(ystar.new,trials_star,p0,log=TRUE))-
    0.5*Det.new-inv.delta*sum(dbinom(ystar.new,trials_star,piML.new,log=TRUE))-
    sum(dbinom(ystar.new,trials_star,p.prop,log=TRUE))

  log.den <- inv.delta*sum(dbinom(ystar.old,trials_star,p.gamma,log=TRUE))+
    1/psi*sum(dbinom(ystar.old,trials_star,p0,log=TRUE))-
    0.5*Det.old-inv.delta*sum(dbinom(ystar.old,trials_star,piML.old,log=TRUE))-
    sum(dbinom(ystar.old,trials_star,p.prop,log=TRUE))

  a.p = min(exp(log.num-log.den),1)

  return (a.p)
}

delta.accept.prob <- function(delta.new,psi.new,delta.old, psi.old,
                              ystar,trials_star,p.gamma,p0,
                              piML_delta ,
                              #                              propA.new,
                              #                              propA.old,
                              propB,include,hyper.type){
  n.star <- length(trials_star)
  p.gamma<-rep(p.gamma,length.out=n.star)
  p0<-rep(p0,length.out=n.star)
  inv.delta.new<-1/delta.new
  inv.delta.old <- 1/delta.old

  if(hyper.type=='gamma')     {
    prior.d.new<-dgamma(delta.new,shape=n.star*0.001,rate=0.001,log=TRUE)
    prior.d.old<-dgamma(delta.old,shape=n.star*0.001,rate=0.001,log=TRUE)}
  if(hyper.type=='ZS')        {
    prior.d.new<- -3/2*log(delta.new)-n/(2*delta.new)
    prior.d.old<- -3/2*log(delta.old)-n/(2*delta.old)}
  if(hyper.type=='hyper-g')   {
    prior.d.new<- -(3/2)*log(1+delta.new)
    prior.d.old<- -(3/2)*log(1+delta.old)
  }
  if(hyper.type=='hyper-g/n') {
    prior.d.new<- -(3/2)*log(1+delta.new/n.star)
    prior.d.old<- -(3/2)*log(1+delta.old/n.star)
  }
  if(hyper.type=='maruyama')  {
    prior.d.new<- ((n-sum(include)+1)/2-7/4)*log(delta.new)-((n-sum(include))/2)*log(1+delta.new)
    prior.d.old<- ((n-sum(include)+1)/2-7/4)*log(delta.old)-((n-sum(include))/2)*log(1+delta.old)
  }

  log.num <-  - 0.5*sum(include)*log(delta.new)+
    inv.delta.new*(sum(dbinom(ystar,trials_star,p.gamma,log=TRUE))-
                     sum(dbinom(ystar,trials_star,piML_delta,log=TRUE)))+
    1/psi.new*sum(dbinom(ystar,trials_star,p0,log=TRUE))-
    #dgamma(delta.new,shape=propA.old,rate=propB,log=TRUE)+
    dlnorm(delta.new,meanlog = log(delta.old),sdlog = propB, log = TRUE)+
    prior.d.new

  log.den <-  - 0.5*sum(include)*log(delta.old)+
    inv.delta.old*(sum(dbinom(ystar,trials_star,p.gamma,log=TRUE))-
                     sum(dbinom(ystar,trials_star,piML_delta,log=TRUE)))+
    1/psi.old*sum(dbinom(ystar,trials_star,p0,log=TRUE))-
    #dgamma(delta.old,shape=propA.new,rate=propB,log=TRUE)+
    dlnorm(delta.old,meanlog = log(delta.new),sdlog = propB, log = TRUE)+
    prior.d.old

  if(log.den>=710){log.den <- 709.78}
  if(log.num>=710){log.num <- 709.78}
  if(log.den<=-710){log.den <- -709.78}
  if(log.num<=-710){log.num <- -709.78}

  a.p = min(exp(log.num-log.den),1)

  return (a.p)
}


pep.glm<-function(y,X,iter,discard,family,n.star=nrow(y), prior.type='CRPEP',
                  model='search',model.prob='beta',hyper=FALSE,
                  hyper.type='hyper-g/n',hyper.multiplier=0.25,ini.b='none',
                  ini.b0='none',ini.y='none',variables=rep(1,ncol(X))){
  Begin<-Sys.time()


  if(model=='single') {
    if(length(variables)!=dim(X)[2]) stop('Length of variables is not equal to the number of columns of design matrix')
    ini.gamma<-c(1,variables)
    include<-ini.gamma==1
    exclude<-ini.gamma==0
  }


  if(family=='binomial') {

    trials<-y[,1]
    success<-y[,2]
    fail<-trials-success
    Y<-cbind(success,fail)
    n<-length(trials)


    if(ini.b0=='none') {
      glm.ini.b0<-glm(Y~1,family=binomial)
      ini.b0<-glm.ini.b0$coefficients
    }

    ini.p0<-expit(ini.b0)

    if(class(ini.y)=='numeric') {
      n.star<-length(ini.y)
      trials.star<-rep(trials,length.out=n.star)
    }
    if(ini.y=='none'){
      trials.star<-rep(trials,length.out=n.star)
      ini.y<-rbinom(n.star,trials.star,ini.p0)
    }
    if(n.star<n) stop('Size of imaginary data must be at least equal to real sample size')

    delta<-n.star
    inv.delta<-1/n.star
    if(prior.type=="CRPEP"){
      psi=1
    }else if(prior.type=="DRPEP"){
      psi=delta
    }

    ini.fail<-trials.star-ini.y
    Y.star<-cbind(ini.y,ini.fail)
    Y.all<-rbind(Y,Y.star)

    X_1<-cbind(1,X)
    if(is.null(colnames(X))) {colnames(X)<-paste(rep('X',ncol(X)),c(1:ncol(X)),sep='')}
    colnames(X_1)<-c('intercept',colnames(X))
    X.star<-X[rep(1:n,length.out=n.star),]
    X.star_1<-cbind(1,X.star)
    X.all<-rbind(X,X.star)

    weights.star<-rep(inv.delta,n.star)
    weights.all<-c(rep(1,n),rep(inv.delta,n.star))

    p<-dim(X_1)[2]

    if(model=='search'){
      ini.gamma<-c(1,rbinom((p-1),1,0.5))
      include<-ini.gamma==1
      exclude<-ini.gamma==0
      act.model<-rep(NA,iter)
      act.predictors<-rep(NA,iter)
    }

    if(model=='single') {
      if(ini.b=='none')   {
        glm.ini.b<-glm(Y.all~X.all[,include[-1]],family=binomial)
        ini.b<-glm.ini.b$coefficients
      }
      ini.b.gamma<-rep(0,p)
      ini.b.gamma[include]<-ini.b
      ini.p.gamma<-expit(X_1%*%ini.b.gamma)
    }


    if(model=='search') {
      if(ini.b=='none')   {
        glm.ini.b<-glm(Y.all~X.all,family=binomial)
        ini.b<-glm.ini.b$coefficients
      }
      ini.b.gamma<-ini.b*ini.gamma
      ini.p.gamma<-expit(X_1%*%ini.b.gamma)
    }

    cand.b<-rep(NA,p)
    cand.b.gamma<-matrix(0,iter,p)
    cand.y<-matrix(NA,iter,n.star)
    if(model=='search'){
      cand.gamma<-matrix(NA,iter+1,p-1)
      cand.gamma[1,]<-ini.gamma[-1]
    }

    if(hyper==TRUE) {
      cand.delta<-rep(NA,iter)
      ratio<-rep(0,4)
      names(ratio)<-c('betas','beta0','y*','delta')
    } else {
      ratio<-rep(0,3)
      names(ratio)<-c('betas','beta0','y*')
    }

    if(model=='single'){
      glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=binomial)
      BetaML<-glm.b.gamma$coefficients
      CovML<-vcov(glm.b.gamma)
    }


    if(model=='search'){
      glm.b.gamma<-glm(Y~X,family=binomial)
      BetaML<-glm.b.gamma$coefficients
      CovML<-vcov(glm.b.gamma)
      CovML.ind<-diag(diag(CovML))
    }

    for(i in 1:iter){
      if (i%%1000 == 0){cat(paste(i," "))}


      #### Generate  betas with gamma ==1 ####


      if(model=='single'){
        cand.b<-rmvnorm(1,BetaML,CovML)
        cand.b.gamma[i,include]<-cand.b
        cand.p.gamma<-expit(X_1%*%cand.b.gamma[i,])

        mh.ratio.b<-b.accept.prob(cand.b.gamma[i, include],cand.p.gamma,
                                  ini.b.gamma[include],ini.p.gamma,success,trials,
                                  ini.y,trials.star,BetaML,CovML,inv.delta)
        #log.posterior.b(cand.b.gamma[i,include],cand.p.gamma,success,fail,ini.y,ini.fail,BetaML,CovML)-
        #log.posterior.b(ini.b.gamma[include],ini.p.gamma,success,fail,ini.y,ini.fail,BetaML,CovML)

        if(runif(1) <= mh.ratio.b) {
          ini.b.gamma<-cand.b.gamma[i,]
          ini.p.gamma<-cand.p.gamma
          ratio[1]<-ratio[1]+1
        } else {
          cand.b.gamma[i,]<-ini.b.gamma
          cand.p.gamma<-ini.p.gamma
          cand.b<-ini.b
        }
      }


      if(model=='search'){

        if(sum(include)==1){
          glm.b.gamma<-glm(Y.all~1,weights=weights.all,family=binomial)
        }
        if(sum(include)!=1){
          glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=binomial)
        }
        BetaML.INC<-glm.b.gamma$coefficients
        CovML.INC<-vcov(glm.b.gamma)
        cand.b<-rmvnorm(1,BetaML.INC,CovML.INC)

        cand.b.gamma[i,include]<-cand.b
        cand.p.gamma<-expit(X_1%*%cand.b.gamma[i,])
        mh.ratio.b<-b.accept.prob(cand.b.gamma[i, include],cand.p.gamma,
                                  ini.b.gamma[include],ini.p.gamma,success,trials,
                                  ini.y,trials.star,BetaML.INC,CovML.INC,inv.delta)
        #log.posterior.b(cand.b,cand.p.gamma,success,fail,ini.y,ini.fail,BetaML.INC,CovML.INC)-
        #log.posterior.b(ini.b.gamma[include],ini.p.gamma,success,fail,ini.y,ini.fail,BetaML.INC,CovML.INC)

        if(runif(1) <= mh.ratio.b) {
          ini.b.gamma<-cand.b.gamma[i,]
          ini.p.gamma<-cand.p.gamma
          ratio[1]<-ratio[1]+1
        } else {
          cand.b.gamma[i,]<-ini.b.gamma
          cand.p.gamma<-ini.p.gamma
          cand.b<-ini.b[include]
        }
      }



      #### Generate betas with gamma ==0 ####

      if(model=='search') {
        cand.b.gamma.exl<-rmvnorm(1,BetaML,CovML.ind)
        cand.b.gamma.exl<-cand.b.gamma.exl[exclude]
      }


      #### Generation and MH for beta0 - NULL model ####
      glm.b0<-glm(Y.star~1,family=binomial)
      BetaML0<-glm.b0$coefficients
      CovML0<-sqrt(psi*vcov(glm.b0))

      cand.b0<-rnorm(1,BetaML0,CovML0)
      cand.p0<-expit(cand.b0)


      mh.ratio.b0<- b0.accept.prob(cand.b0,ini.b0,ini.y,trials.star,BetaML0,CovML0,psi)

      if(runif(1) <= mh.ratio.b0)  {
        ini.b0<-cand.b0
        ini.p0<-cand.p0
        ratio[2]<-ratio[2]+1
      } else {
        cand.b0<-ini.b0
        cand.p0<-ini.p0
      }


      #### Generation and MH for y* ####
      logy <- 1/psi*log(1-cand.p0)+inv.delta*(1-cand.p.gamma)
      logx <- 1/psi*log(cand.p0)+inv.delta*(cand.p.gamma)

      p.STAR<- expit(logx-logy)
      p.STAR<-rep(p.STAR,length.out=n.star)

      cand.y[i,]<-rbinom(n.star,trials.star,p.STAR)
      cand.fail<-trials.star-cand.y[i,]
      Y.star.cand<-cbind(cand.y[i,],cand.fail)

      # Laplace approximation at t-1
      if(sum(include)==1){
        glm.prev<-glm(Y.star~1,family=binomial)
      }else if(sum(include)!=1){
        glm.prev<-glm(Y.star~X.star[,include[-1]],family=binomial)
      }
      betaML.prev<-glm.prev$coefficients
      piML.prev <- glm.prev$fitted.values
      H.prev<-diag(c(trials.star*piML.prev*(1-piML.prev)))
      Det.prev<--as.numeric(determinant(inv.delta*(t(X.star_1[,include])%*%H.prev%*%X.star_1[,include]),log=TRUE)$modulus)


      # Laplace approximation at t
      if(sum(include)==1){
        glm.curr<-glm(Y.star.cand~1,family=binomial)
      }else if(sum(include)!=1){
        glm.curr<-glm(Y.star.cand~X.star[,include[-1]],family=binomial)
      }
      betaML.curr<-glm.curr$coefficients
      piML.curr <- glm.curr$fitted.values
      H.curr<-diag(c(trials.star*piML.curr*(1-piML.curr)))
      Det.curr<--as.numeric(determinant(inv.delta*(t(X.star_1[,include])%*%H.curr%*%X.star_1[,include]),log=TRUE)$modulus)


      mh.ratio.y<- ystar.accept.prob(cand.y[i,],ini.y, trials.star,
                                     cand.p.gamma, cand.p0, piML.curr,
                                     piML.prev, Det.curr, Det.prev, p.STAR,
                                     psi,inv.delta)

      if(runif(1) <= mh.ratio.y) {
        ini.y<-cand.y[i,]
        ini.fail<-cand.fail
        Y.star<-Y.star.cand
        ratio[3]<-ratio[3]+1
      } else {
        cand.y[i,]<-ini.y
        cand.fail<-ini.fail
        Y.star.cand<-Y.star
      }

      Y.all<-rbind(Y,Y.star)

      #### Generation and MH for delta ####
      if(hyper==TRUE){
        #A.prev<-delta*hyper.multiplier
        B<-hyper.multiplier
        cand.delta[i]<-  rlnorm(1,meanlog = log(delta),sdlog = B)#rgamma(1,shape=A.prev,rate=B)#
        cand.inv.delta<-1/cand.delta[i]
        #A.curr<-cand.delta[i]*hyper.multiplier
        if(prior.type=="CRPEP"){
          cand.psi <- 1
        }else if(prior.type=="DRPEP"){
          cand.psi <- cand.delta[i]
        }
        if(sum(include)==1){
          glm.delta<-glm(Y.star~1,family=binomial)
        }
        if(sum(include)!=1){
          glm.delta<-glm(Y.star~X.star[,include[-1]],family=binomial)
        }
        betaML.delta<-glm.delta$coefficients
        piML.delta<-glm.delta$fitted.values


        mh.ratio.delta<- delta.accept.prob(cand.delta[i],cand.psi,delta, psi,
                                           ini.y,trials.star,cand.p.gamma,cand.p0,
                                           piML.delta ,
                                           #A.curr, A.prev,
                                           B,include,hyper.type)

        if(runif(1) <= mh.ratio.delta) {
          inv.delta<-cand.inv.delta
          delta<-cand.delta[i]
          if(prior.type=="CRPEP"){psi <- 1}else if(prior.type=="DRPEP"){psi <- cand.psi}
          ratio[4]<-ratio[4]+1
        } else {
          cand.inv.delta<-inv.delta
          cand.delta[i]<-delta
          if(prior.type=="CRPEP"){psi <- 1}else if(prior.type=="DRPEP"){cand.psi <- psi}
        }

        #weights.star<-rep(inv.delta,n.star)
        #weights.all<-c(rep(1,n),rep(inv.delta,n.star))
        if(delta<0.001){print(paste("delta=",delta))}

      }


      #### Generation of gamma's ####
      if(model=="search"){

        cand.b.all<-cand.b.gamma[i,]
        cand.b.all[exclude]<-cand.b.gamma.exl

        for(j in 2:p) {

          inc.state<-rbind(include,include)

          inc.state[1,j]<-1
          inc.state[2,j]<-0

          exc.state<-inc.state
          exc.state<-inc.state+1
          exc.state[exc.state==2]<-0

          cand.b.inc1<-c(cand.b.all*inc.state[1,])
          cand.b.inc0<-c(cand.b.all*inc.state[2,])

          cand.b.exc1<-c(cand.b.all*exc.state[1,])
          cand.b.exc0<-c(cand.b.all*exc.state[2,])

          Pi_y1<-expit(X_1%*%cand.b.inc1)
          Pi_y0<-expit(X_1%*%cand.b.inc0)
          Lik_y1<-sum(dbinom(success,trials,Pi_y1,log=TRUE))
          Lik_y0<-sum(dbinom(success,trials,Pi_y0,log=TRUE))

          Lik_y.star1<-inv.delta*sum(dbinom(cand.y[i,],trials.star,Pi_y1,log=TRUE))
          Lik_y.star0<-inv.delta*sum(dbinom(cand.y[i,],trials.star,Pi_y0,log=TRUE))

          inc.ind1<-inc.state[1,]==1
          inc.ind0<-inc.state[2,]==1

          glm.inc1<-glm(Y.star.cand~X.star[,inc.ind1[-1]],family=binomial)
          betaML.inc1<-glm.inc1$coefficients
          piML.inc1<-glm.inc1$fitted.values
          H.inc1<-diag(c(trials.star*piML.inc1*(1-piML.inc1)))
          Det.inc1<--as.numeric(determinant(inv.delta*(t(X.star_1[,inc.ind1])%*%H.inc1%*%X.star_1[,inc.ind1]),logarithm = TRUE)$modulus)

          if(sum(inc.ind0)==1){
            glm.inc0<-glm(Y.star.cand~1,family=binomial)
          }
          if(sum(inc.ind0)!=1){
            glm.inc0<-glm(Y.star.cand~X.star[,inc.ind0[-1]],family=binomial)
          }
          betaML.inc0<-glm.inc0$coefficients
          piML.inc0<-glm.inc0$fitted.values
          H.inc0<-diag(c(trials.star*piML.inc0*(1-piML.inc0)))
          Det.inc0<--as.numeric(determinant(inv.delta*(t(X.star_1[,inc.ind0])%*%H.inc0%*%X.star_1[,inc.ind0]),logarithm = TRUE)$modulus)

          Marg1<-sum(inc.state[1,])/2*log(2*pi)+0.5*Det.inc1+inv.delta*sum(dbinom(cand.y[i,],trials.star,piML.inc1,log=TRUE))
          Marg0<-sum(inc.state[2,])/2*log(2*pi)+0.5*Det.inc0+inv.delta*sum(dbinom(cand.y[i,],trials.star,piML.inc0,log=TRUE))

          exc.ind1<-exc.state[1,]==1
          exc.ind0<-exc.state[2,]==1


          if(sum(exc.ind1)==0){Pseudo1<-0}
          if(sum(exc.ind1)==1){Pseudo1<-dnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],sqrt(CovML.ind[exc.ind1,exc.ind1]),log=TRUE)}
          if(sum(exc.ind1)>1){Pseudo1<-dmvnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],CovML.ind[exc.ind1,exc.ind1],log=TRUE)}

          if(sum(exc.ind0)==1){Pseudo0<-dnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],sqrt(CovML.ind[exc.ind0,exc.ind0]),log=TRUE)}
          if(sum(exc.ind0)>1){Pseudo0<-dmvnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],CovML.ind[exc.ind0,exc.ind0],log=TRUE)}

          d<-sum(inc.state[1,-1])-1
          if(model.prob=='uniform') {logOj<-Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0}
          if(model.prob=='beta') {logOj<-Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0+log(d+1)-log(p-1-d)}
          prop.gamma.j<-expit(logOj)

          ini.gamma[j]<-rbinom(1,1,prop.gamma.j)
          include<-ini.gamma==1
          exclude<-ini.gamma==0

        } # end of inner for loop for the gamma's

        cand.gamma[i+1,]<-ini.gamma[2:p]

        ini.b.gamma<-cand.b.all*ini.gamma
        ini.p.gamma<-expit(X_1%*%ini.b.gamma)

        act.model[i]<-sum(cand.gamma[i+1,]*2^(0:(p-2)))
        act.predictors[i]<-c(paste0("X", c(1:(p-1))[include[-1]],collapse='+'))
      }# end of if model=='search'
    }# end of iterations


  }# end of binomial

  End<-Sys.time()
  runtime<-difftime(End,Begin)
  if(model=='single'){
    final.b<-cand.b.gamma[(discard+1):iter,include]
    colnames(final.b)<-colnames(X_1[,include])
    if(hyper==FALSE) {results<-list(betas=final.b,acceptance.ratios=ratio/iter,runtime=runtime)}
    if(hyper==TRUE) {
      final.delta<-cand.delta[(discard+1):iter]
      results<-list(betas=final.b,delta=final.delta,acceptance.ratios=ratio/iter,runtime=runtime)
    }
  }
  if(model=='search'){
    final.b<-cand.b.gamma[(discard+1):iter,]
    colnames(final.b)<-colnames(X_1)
    final.gamma<-cand.gamma[(discard+1):iter,]
    colnames(final.gamma)<-colnames(X)
    PIPs<-round(apply(final.gamma,2,sum)/(iter-discard),3)
    act.model<-act.model[(discard+1):iter]
    act.predictors<-act.predictors[(discard+1):iter]
    crosstab<-table(act.model)
    rows<-c()
    for(i in 1:length(crosstab)){
      rows[i]<-which(act.model==names(crosstab)[i])[1]
    }
    summary<-data.frame(crosstab,act.predictors[rows])
    summary<-summary[order(summary[,2],decreasing=TRUE),]
    rownames(summary)<-1:dim(summary)[1]
    summary<-data.frame(summary[,1:2],summary[,2]/(iter-discard),summary[,3])
    names(summary)<-c('Model','Frequency','PostProb','Variables')
    if(hyper==FALSE) {results<-list(betas=final.b,gammas=final.gamma,search=summary,PIPs=PIPs,acceptance.ratios=ratio/iter,runtime=runtime)}
    if(hyper==TRUE) {
      final.delta<-cand.delta[(discard+1):iter]
      results<-list(betas=final.b,gammas=final.gamma,delta=final.delta,search=summary,PIPs=PIPs,acceptance.ratios=ratio/iter,runtime=runtime)
    }
  }
  results

}
