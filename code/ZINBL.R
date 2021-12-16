################
#ZINBLDA
#estimated parameter based on the ZINBLDA model

Soft <- function(x,a){
  return(sign(x)*pmax(abs(x)-a,0))
}

fun <- function(x,mu,disperhatmean){
  signx <- sign(x==0)
  zinb <- function(p) {
    
    res <- sum(sign(x==0)*log(p[1]+(1-p[1])*(1/(1+p[3]*p[2]))^(1/p[2]))+
                 (1-sign(x==0))*(log(1-p[1])+lgamma(x+1/p[2])-lgamma(x+1)-
                                   lgamma(1/p[2])+x*log(p[3]*p[2])-x*log(1+p[3]*p[2])-(1/p[2])*log(1+p[3]*p[2])))
    
    return(-res)
    
  }
  nlminb(c(0,1,mu),zinb,lower=c(0.01,0.001,0.01),upper=c(0.999,disperhatmean+0.5,9999),control = list(step.max=0.2))$par
}



permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}


GetD <- function(ns, x, y, rho,beta,rhos=NULL){
  if(!is.null(rho) && !is.null(rhos)) stop("do you want to use rho or rhos in GetD function???")
  if(is.null(rhos)){
    uniq <- sort(unique(y))
    ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
    for(k in 1:length(uniq)){
      a <- colSums(x[y==uniq[k],])+beta
      b <- colSums(ns[y==uniq[k],])+beta
      ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
    }
    return(ds)
  } else {
    uniq <- sort(unique(y))
    ds.list <- list()
    for(rho in rhos){
      ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        ds[k,] <- 1+Soft(a/b-1,rho/sqrt(b))
      }
      ds.list[[which(rhos==rho)]] <- ds
    }
    return(ds.list)
  }
}


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}




ZINB.cv <-
  function(x,y,rhos=NULL,beta=1,nfolds=5,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"),folds=NULL, prior=NULL){
    type <- match.arg(type)
    if(is.null(rhos)){

      ns <- NullModel(x,type=type)$n
      uniq <- sort(unique(y))
      maxrho <- rep(NA, length(uniq))
      for(k in 1:length(uniq)){
        a <- colSums(x[y==uniq[k],])+beta
        b <- colSums(ns[y==uniq[k],])+beta
        maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
      }
      rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
    }
    if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
    nfolds <- length(folds)
    errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
    for(i in 1:nfolds){
      cat(i,fill=FALSE)
      tr <- -folds[[i]]
      te <- folds[[i]]
      out <-ZINBLDA(x[tr,],y[tr],x[te,],rhos=rhos,phihat=phihat,beta=beta,prob0=prob0,type="mle", prior=prior) # Have already power-transformed x, so don't need to do it again!!!
      for(j in 1:length(rhos)){      
        errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
        nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
      }
    }
    cat(fill=TRUE)
    save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds,type=type)
    return(save)
  }


ZINBLDA<-
  function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,phihat=0,prob0=NULL,type=c("mle","deseq","quantile"), prior=NULL){
    eps=1e-8
    if(is.null(xte)){
      xte <- x
      warning("Since no xte was provided, testing was performed on training data set.")
    }
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx3<-sign(xte==0)
    
    if(is.null(rhos)){
      ds <- GetD(ns,x,y,rho,beta)
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        for(i in 1:nrow(xte)){
          
          dstar = ds[k,]
          part2=nste[i,]*dstar 
          part3<-(1/(1+part2*phihat))^((1+eps)/(phihat+eps))
          discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/(phihat+eps))*log(1+part2*phihat))+log(prior[k]+1)
        }
      }
      save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type)
      return(save)
    } else {
      save <- list()
      ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
      for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in 1:length(uniq)){
          for(i in 1:nrow(xte))   {
            
            dstar = ds[k,]
            part2=nste[i,]*dstar 
            part3<-(1/(1+part2*phihat))^(1/phihat)
            discriminant[i,k] <-sum(signx3[i,]*log(prob0[i,]+(1-prob0[i,])*part3+1))+sum(xte[i,]*(1-signx3[i,])*(log(dstar)-log(1+part2*phihat)))-sum((1-signx3[i,])*(1/(phihat+eps))*log(1+part2*phihat))+log(prior[k])
          }
        }
        save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=apply(discriminant,1,which.max),rho=rho,x=x,y=y,xte=xte,type=type))
      }
      return(save)
    }
  }

select2<-function(dat,gene_no_list){
  n=length(dat$yte)
  K=2
  dddnum <- 1
  sdddnum<- 1
  c1=sum(dat$y==1)
  c2=sum(dat$y==2)
  cdata<-t(dat$x)
  ctedata<-t(dat$xte)
  data<-cbind(t(dat$x),t(dat$xte))
  fW <- calcNormFactors(cdata)
  for(l in 1:length(fW)){
    if(is.na(fW[l]))
      fW[l]<-1
  }
  newdata<-NULL
  for(i in 1:length(fW)){
    newdata<-cbind(newdata,cdata[,i]/fW[i])
  }
  a1<-newdata[,dat$y==1]
  a2<-newdata[,dat$y==2]
  trainx=cbind(a1,a2)
  
  myscore <- rep(0,length=nrow(trainx))####
  smyscore <- rep(0,length=nrow(trainx))
  newylabel<-matrix(0,nrow=nrow(trainx),ncol = n)
  classmean<-1:K
  
  
  for (i in 1:nrow(trainx)){ 
    for(iii in 1:K){
      classmean[iii]=mean(newdata[i,dat$y==iii])
    }
    crank= rank(classmean)
    newylabel[i,]<- c(rep(crank[1],c1),rep(crank[2],c2))
    rank<-rank(trainx[i,])
    rankclass<-rank(newylabel[i,])
    sdIX=sd(rank)
    sdIY=sd(rankclass)
    sum<-0
    for(jj in 1:n){
      sum=sum+(rank[jj]-rankclass[jj])^2
    }
    fenzi<-(sdIX^2+sdIY^2-sum/n)/2
    smyscore[i]=fenzi/(sdIX*sdIY)
    rmean <- mean(newdata[i,])
    BSS<-0
    WSS <-0.01
    for(jj in 1:K){
      rc1mean <- mean(newdata[i,dat$y==jj])
      BSS <-BSS+ sum(dat$y==jj)*(rc1mean - rmean)^2 
      WSS <- WSS + sum((newdata[i,dat$y==jj] -rmean)^2)
    }
    myscore[i] <- BSS/WSS
    
    
  }
  sort(smyscore,decreasing=TRUE)
  sorttrainx <- sort.list(myscore, decreasing=TRUE)
  sorttrainx <- sorttrainx[1:gene_no_list]
  ddd=sorttrainx
  data1<-data[ddd,]
  
  
  x=t(data1[,1:n])
  xte=t(data1[,(n+1):(2*n)])
  dat <- list(x=x,y=dat$y,xte=xte,yte=dat$yte)
  return (dat)
}


select3<-function(dat,gene_no_list){
  
  n=length(dat$yte)
  K=3
  dddnum <- 1
  sdddnum<- 1
  c1=sum(dat$y==1)
  c2=sum(dat$y==2)
  c3=sum(dat$y==3)
  cdata<-t(dat$x)
  ctedata<-t(dat$xte)
  data<-cbind(t(dat$x),t(dat$xte))
  fW <- calcNormFactors(cdata)
  for(l in 1:length(fW)){
    if(is.na(fW[l]))
      fW[l]<-1
  }
  newdata<-NULL
  for(i in 1:length(fW)){
    newdata<-cbind(newdata,cdata[,i]/fW[i])
  }
  a1<-newdata[,dat$y==1]
  a2<-newdata[,dat$y==2]
  a3<-newdata[,dat$y==3]
  trainx=cbind(a1,a2,a3)
  
  myscore <- rep(0,length=nrow(trainx))####
  smyscore <- rep(0,length=nrow(trainx))
  newylabel<-matrix(0,nrow=nrow(trainx),ncol = n)
  classmean<-1:K
  
  
  for (i in 1:nrow(trainx)){ 
    for(iii in 1:K){
      classmean[iii]=mean(newdata[i,dat$y==iii])
    }
    crank= rank(classmean)
    newylabel[i,]<- c(rep(crank[1],c1),rep(crank[2],c2),rep(crank[3],c3))
    rank<-rank(trainx[i,])
    rankclass<-rank(newylabel[i,])
    sdIX=sd(rank)
    sdIY=sd(rankclass)
    sum<-0
    for(jj in 1:n){
      sum=sum+(rank[jj]-rankclass[jj])^2
    }
    fenzi<-(sdIX^2+sdIY^2-sum/n)/2
    smyscore[i]=fenzi/(sdIX*sdIY)
    rmean <- mean(newdata[i,])
    BSS<-0
    WSS <-0.01
    for(jj in 1:K){
      rc1mean <- mean(newdata[i,dat$y==jj])
      BSS <-BSS+ sum(dat$y==jj)*(rc1mean - rmean)^2 
      WSS <- WSS + sum((newdata[i,dat$y==jj] -rmean)^2)
    }
    myscore[i] <- BSS/WSS
    
    
  }
  sort(smyscore,decreasing=TRUE)
  sorttrainx <- sort.list(myscore, decreasing=TRUE)
  sorttrainx <- sorttrainx[1:gene_no_list]
  ddd=sorttrainx
  data1<-data[ddd,]
  
  
  x=t(data1[,1:n])
  xte=t(data1[,(n+1):(2*n)])
  dat <- list(x=x,y=dat$y,xte=xte,yte=dat$yte)
  return (dat)
}

select4<-function(dat,gene_no_list){
  
  n=length(dat$yte)
  K=4
  dddnum <- 1
  sdddnum<- 1
  c1=sum(dat$y==1)
  c2=sum(dat$y==2)
  c3=sum(dat$y==3)
  c4=sum(dat$y==4)
  cdata<-t(dat$x)
  ctedata<-t(dat$xte)
  data<-cbind(t(dat$x),t(dat$xte))
  fW <- calcNormFactors(cdata)
  for(l in 1:length(fW)){
    if(is.na(fW[l]))
      fW[l]<-1
  }
  newdata<-NULL
  for(i in 1:length(fW)){
    newdata<-cbind(newdata,cdata[,i]/fW[i])
  }
  a1<-newdata[,dat$y==1]
  a2<-newdata[,dat$y==2]
  a3<-newdata[,dat$y==3]
  a4<-newdata[,dat$y==4]
  trainx=cbind(a1,a2,a3,a4)
  
  myscore <- rep(0,length=nrow(trainx))####
  smyscore <- rep(0,length=nrow(trainx))
  newylabel<-matrix(0,nrow=nrow(trainx),ncol = n)
  classmean<-1:K
  
  
  for (i in 1:nrow(trainx)){ 
    for(iii in 1:K){
      classmean[iii]=mean(newdata[i,dat$y==iii])
    }
    crank= rank(classmean)
    newylabel[i,]<- c(rep(crank[1],c1),rep(crank[2],c2),rep(crank[3],c3),rep(crank[4],c4))
    rank<-rank(trainx[i,])
    rankclass<-rank(newylabel[i,])
    sdIX=sd(rank)
    sdIY=sd(rankclass)
    sum<-0
    for(jj in 1:n){
      sum=sum+(rank[jj]-rankclass[jj])^2
    }
    fenzi<-(sdIX^2+sdIY^2-sum/n)/2
    smyscore[i]=fenzi/(sdIX*sdIY)
    rmean <- mean(newdata[i,])
    BSS<-0
    WSS <-0.01
    for(jj in 1:K){
      rc1mean <- mean(newdata[i,dat$y==jj])
      BSS <-BSS+ sum(dat$y==jj)*(rc1mean - rmean)^2 
      WSS <- WSS + sum((newdata[i,dat$y==jj] -rmean)^2)
    }
    myscore[i] <- BSS/WSS
    
    
  }
  sort(smyscore,decreasing=TRUE)
  sorttrainx <- sort.list(myscore, decreasing=TRUE)
  sorttrainx <- sorttrainx[1:gene_no_list]
  ddd=sorttrainx
  data1<-data[ddd,]
  
  
  x=t(data1[,1:n])
  xte=t(data1[,(n+1):(2*n)])
  dat <- list(x=x,y=dat$y,xte=xte,yte=dat$yte)
  return (dat)
}

ZINBScore<-function(dat,batch_size){
  m<-DeepZINB(X_train_=dat$x,y_train=dat$y,X_test_=dat$xte,y_test=dat$yte,batch_size=batch_size)
  y_tr<-array(unlist(y_train))
  y_te<-array(unlist(y_test))
  n_disperhatZINB=colMeans(m$disp)
  n_disperhatZINB_t=colMeans(m$disp_t)
  n_prob0_ZINB_train<-m$pi
  n_prob0_ZINB_test<-m$pi_t
  DeepZINBLcv.out <- ZINB.cv(x=m$X_train,y=y_tr,prob0=n_prob0_ZINB_train,
                             phihat=n_disperhatZINB)
  DeepZInbL.out <- ZINBLDA(x=m$X_train,y=y_tr,xte=m$X_test,phihat=n_disperhatZINB_t,
                           prob0=n_prob0_ZINB_test,rho=DeepZINBLcv.out$bestrho) 
  y_pre=DeepZInbL.out$ytehat
  return(list(y_pre=y_pre,y_test=y_test))
}

