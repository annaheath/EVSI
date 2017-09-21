##evsi.calc############################################################
evsi.calc<-function(comp.evsi.N,wtp=NULL,N=NULL,CI=NULL){
  ##'Calculates the EVSI accross willingness-to-pay, sample size and gives uncertainty for these
  ##'measurements as required. This function can be fed into/used in the other plot functions
  ##'This should allow people to upload their own EVSI matrices and still use the plotting capabilities.
  ##INPUTS
  ##'@param comp.evsi.N Output of the comp.evsi.N function
  ##'@param wtp The willingness to pay value(s) for which the EVSI should be calculated.
  ##'   It will be chosen as all the wtp values if NULL.
  ##'@param N The sample sizes for which the EVSI should be calculated. If NULL it will select
  ##'   all the values that the EVSI has been calculated for. NOTE: N can be chosen as any
  ##'   value but wtp needs to be on the grid chosen in comp.evsi.N.
  ##'@param CI The confidence levels for the credible intervals of the EVSI. If NULL chosen as
  ##'
  ##'OUTPUTS
  ##'@return An evsi object.
  ##'1. evsi An array containing the EVSI by wtp, N and across different uncertaincies
  ##'2. attrib A list of wtp, N and prob describing the attributes of the evsi matrix.
  ##'3. evppi An evppi object containing all the information about the calculation of
  ##'   EVPPI.
  ##'4. he A bcea object containing all the information about the underlying health
  ##'   economic model
  ##'   @example
  ##'   ...

  if(class(comp.evsi.N)!="evsi.N"){stop("comp.evsi.N must be calculated using the comp.evsi.N function. Please use this function. evsi objects can be created using the evsi.upload function.")}
  #Format wtp
  cl<-class(wtp)
  #Select all wtp values if NULL
  if(cl!="numeric"){
    wtp<-comp.evsi.N$wtp
  }
  wtp.length<-length(wtp)
  #Select wtp from the grid of wtp if chosen by user.
  if(cl=="numeric"){
    for(i in 1:wtp.length){
      wtp[i]<-comp.evsi.N$wtp[which.min(abs(wtp[i]-comp.evsi.N$wtp))]
    }
  }

  #Format N
  cl<-class(N)
  #Select all N values if NULL
  if(is.null(N)){
    N<-seq(min(comp.evsi.N$N),max(comp.evsi.N$N),by=10)
  }
  N.length<-length(N)

  #Format CI
  cl<-class(CI)
  #Set the default CI values if NULL
  if(cl!="numeric"){
    CI<-c(0.025,0.25,0.5,0.75,0.975)
  }
  CI.length<-length(CI)
  #Set N as requested by user
  if(cl=="numeric"){
    CI<-CI
  }
  #Use CI to create a beta matrix
  #Find the appropriate quantiles for the beta parameter of interest
  beta.quantiles<-apply(comp.evsi.N$beta,c(2,3),quantile,prob=CI)

  #Extracting the beta parameter for the wtp of interest
  beta.focal<-beta.quantiles[,which(comp.evsi.N$wtp%in%wtp),]

  if(length(dim(beta.focal))==2){dim(beta.focal)<-c(dim(beta.focal),1)}
  EVSI.mat<-array(NA,dim=c(N.length,wtp.length,CI.length))
  INB.full<-array(comp.evsi.N$he$ib,c(length(comp.evsi.N$he$k),comp.evsi.N$he$n.sim,comp.evsi.N$he$n.comparisons))

  calc.EVSI<-function(beta.focal,wtp,N){
    INB.full<--(wtp*comp.evsi.N$he$delta.e-comp.evsi.N$he$delta.c)
    var.full<-as.matrix(var(INB.full))
    #Find the fitted values for the wtp
    INB<--(wtp*comp.evsi.N$evi$fitted.effects[,-comp.evsi.N$he$n.comparators]-
      comp.evsi.N$evi$fitted.costs[,-comp.evsi.N$he$n.comparators])
    var.INB<-var(INB)
    #Variance for a specific N and beta
    var.pre<-var(INB)*(N/(N+beta.focal))

    if(comp.evsi.N$he$n.comparisons==1){
      pre.post.var<-max(0,var.pre)
      INB.star<-(INB-mean(INB))/sd(INB)*sqrt(pre.post.var)+mean(INB)
      }
    #Rescaled fitted values
    if(comp.evsi.N$he$n.comparisons>1){pre.post.var<-var.pre
    check<-base::eigen(pre.post.var)
    #Defintion SQRT from https://stat.ethz.ch/pipermail/r-help/2007-January/124147.html
    pre.post.var.sqrt<-check$vectors%*%diag(sqrt(pmax(0,check$values)))%*%t(check$vectors)
    INB.mean<-matrix(rep(base::colMeans(INB),comp.evsi.N$he$n.sim),nrow=comp.evsi.N$he$n.sim,byrow=TRUE)

    #Fast matrix square root inverse
    decom<-eigen(var.INB)
    var.INB.sqrt.inv<-chol2inv(chol(decom$vectors%*%diag(sqrt(decom$values))%*%t(decom$vectors)))

    #Rescale fitted INB
    INB.star<-  (INB-INB.mean)%*%var.INB.sqrt.inv%*%
                       pre.post.var.sqrt+INB.mean}
    #Calculate EVSI
    EVSI<-mean(apply(cbind(INB.star,0),1,max))-max(apply(cbind(INB.star,0),2,mean))
    return(EVSI)
  }

  n.choices<-1
  index<-matrix(NA,nrow=comp.evsi.N$he$n.comparisons,ncol=comp.evsi.N$he$n.comparisons)
  for(i in 1:comp.evsi.N$he$n.comparisons){
    for(j in i:comp.evsi.N$he$n.comparisons){
      index[i,j]<-index[j,i]<-n.choices
      n.choices<-n.choices+1
    }
  }

  for(i in 1:CI.length){
    for(j in 1:wtp.length){
      beta.of.interest<-matrix(beta.focal[i,j,index],nrow=comp.evsi.N$he$n.comparisons,ncol=comp.evsi.N$he$n.comparisons)
      EVSI.mat[,j,i]<-sapply(N,calc.EVSI,beta.focal=beta.of.interest,wtp=wtp[j])
    }
  }
  to.return<-list(evsi=EVSI.mat,
                  attrib=list(wtp=wtp,N=N,CI=CI),
                  evppi=comp.evsi.N$evi,
                  he=comp.evsi.N$he)
  class(to.return)<-"evsi"
  return(to.return)
}
