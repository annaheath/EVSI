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
    wtp<-comp.evsi.N$he$k
  }
  wtp.length<-length(wtp)

  #Format N
  cl<-class(N)
  #Select all N values if NULL
  if(is.null(N)){
    N.min<-min(comp.evsi.N$N)
    N.max<-max(comp.evsi.N$N)
    if(N.min==N.max){
      N<-N.min
    }

    N<-unique(round(seq(N.min,N.max,length.out=100)))
  }
  N.length<-length(N)

  #Format CI
  cl<-class(CI)
  #Set the default CI values if NULL
  if(cl!="numeric"){
    CI<-c(0.025,0.25,0.5,0.75,0.975)
  }
  #Set N as requested by user
  if(cl=="numeric"){
    CI<-CI
    CI.inv<-1-CI
    CI<-ordered(unique(CI,CI.inv))
  }
  CI.length<-length(CI)
  #Use CI to create a beta matrix
  #Find the appropriate quantiles for the beta parameter of interest
  beta.quantiles<-array(apply(comp.evsi.N$beta,c(2,3),quantile,prob=CI),dim=c(CI.length,2,dim(comp.evsi.N$beta)[3]))

  #Extracting the beta parameter for the costs and effects
  beta.focal.e<-beta.quantiles[,1,]
  beta.focal.c<-beta.quantiles[,2,]

  e.full<-(comp.evsi.N$he$delta.e)
  c.full<-(comp.evsi.N$he$delta.c)
  e.fit<-comp.evsi.N$evppi$fitted.effects[,-comp.evsi.N$he$n.comparators]
  c.fit<-comp.evsi.N$evppi$fitted.costs[,-comp.evsi.N$he$n.comparators]

  #How to populate the variance matrices
  n.choices<-1
  index<-matrix(NA,nrow=comp.evsi.N$he$n.comparisons,ncol=comp.evsi.N$he$n.comparisons)
  for(i in 1:comp.evsi.N$he$n.comparisons){
    for(j in i:comp.evsi.N$he$n.comparisons){
      index[i,j]<-index[j,i]<-n.choices
      n.choices<-n.choices+1
    }
  }

  #Function to find rescaled costs and effects
  calc.fit<-function(beta.focal,N,fit,full){
    var.full<-as.matrix(var(full))
    #Find the fitted values for the wtp
    var.fit<-var(fit)
    #Variance for a specific N and beta
    var.pre<-var.fit*(N/(N+beta.focal))

    if(comp.evsi.N$he$n.comparisons==1){
      pre.post.var<-max(0,var.pre)
      fit.star<-(fit-mean(fit))/sd(fit)*sqrt(pre.post.var)+mean(fit)
    }
    #Rescaled fitted values
    if(comp.evsi.N$he$n.comparisons>1){
      pre.post.var<-var.pre
      check<-base::eigen(pre.post.var)
      #Defintion SQRT from https://stat.ethz.ch/pipermail/r-help/2007-January/124147.html
      pre.post.var.sqrt<-check$vectors%*%diag(sqrt(pmax(0,check$values)))%*%t(check$vectors)
      fit.mean<-matrix(rep(base::colMeans(fit),comp.evsi.N$he$n.sim),nrow=comp.evsi.N$he$n.sim,byrow=TRUE)

      #Fast matrix square root inverse
      decom<-eigen(var.fit)
      var.full.sqrt.inv<-chol2inv(chol(decom$vectors%*%diag(sqrt(decom$values))%*%t(decom$vectors)))

      #Rescale fitted INB
      fit.star<-  (fit-fit.mean)%*%var.full.sqrt.inv%*%
        pre.post.var.sqrt+fit.mean}
    return(fit.star)
  }

  e.star<-list()
  for(i in 1:CI.length){
    e.star[[i]]<-lapply(N,calc.fit,beta.focal=matrix(as.matrix(beta.focal.e)[i,index],
                                                     nrow=comp.evsi.N$he$n.comparisons,ncol=comp.evsi.N$he$n.comparisons),
                        fit=e.fit,full=e.full)

  }

  c.star<-list()
  for(i in 1:CI.length){
    c.star[[i]]<-lapply(N,calc.fit,beta.focal=matrix(as.matrix(beta.focal.c)[i,index],
                                                     nrow=comp.evsi.N$he$n.comparisons,ncol=comp.evsi.N$he$n.comparisons),
                        fit=c.fit,full=c.full)

  }

  wtp.func<-function(wtp.s){
    INB.star<-wtp.s*(-e.star[[i]][[j]])+c.star[[i]][[j]]
    EVSI<-sum(do.call(pmax,as.data.frame(cbind(INB.star,0))))/comp.evsi.N$he$n.sim-
      max(apply(cbind(INB.star,0),2,function(x){sum(x)/comp.evsi.N$he$n.sim}))
    return(EVSI)
  }

  EVSI.mat<-array(NA,dim=c(N.length,wtp.length,CI.length))
  for(i in 1:CI.length){
    for(j in 1:N.length){
        EVSI.mat[j,,i]<-sapply(wtp,wtp.func)
    }
  }

  to.return<-list(evsi=EVSI.mat,
                  attrib=list(wtp=wtp,N=N,CI=CI),
                  evppi=comp.evsi.N$evppi,
                  he=comp.evsi.N$he)
  class(to.return)<-"evsi"
  return(to.return)
}
