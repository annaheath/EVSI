##evsi.calc############################################################
evsi.calc<-function(mm.var, wtp=NULL, N=NULL, CI=NULL){
  ##'Calculates the EVSI accross willingness-to-pay, sample size and gives uncertainty for these
  ##'measurements as required. This function can be fed into/used in the other plot functions
  ##'This should allow people to upload their own EVSI matrices and still use the plotting capabilities.
  ##INPUTS
  ##'@param mm.var An mm.var object computed from the mm.post.var function in the EVSI package.
  ##'   The posterior variances can also be calculated external to the package and passed as a list
  ##'   to the evsi.calc function. If this is the case a bcea object, an evppi object, N.calc and
  ##'   a Bayesian updating program, if mm.var is estimated across different sample sizes must 
  ##'   also be provided.
  ##'@param wtp The willingness to pay value(s) for which the EVSI should be calculated.
  ##'   It will be chosen as all the wtp values in the bcea object, either passed as he or 
  ##'   directly from the mm.var object.
  ##'@param N The sample sizes for which the EVSI should be calculated. If NULL it will select
  ##'   a range of values for which the EVSI was calculated in the mm.var object. If the mm.var
  ##'   object is not given then this argument must be provided.
  ##'@param CI If EVSI is calculated across multiple sample sizes then this gives the confidence 
  ##'   level for the credible intervals of the EVSI. If NULL chosen as 95% and 75% intervals.
  ##'
  ##OUTPUTS
  ##'@return An evsi object.
  ##'1. evsi An array containing the EVSI by wtp, N and across different uncertaincies
  ##'2. attrib A list of wtp, N and prob describing the attributes of the evsi matrix.
  ##'3. evppi An evppi object containing all the information about the calculation of
  ##'   EVPPI.
  ##'4. he A bcea object containing all the information about the underlying health
  ##'   economic model
  ##'   @example
  ##'   ...


  if(class(mm.var)!= c("mm.var")){stop("mm.var must either be calculated using the mm.post.var function.")}
  
  ### Set up
  # Format wtp
  cl <- class(wtp)
  #Select all wtp values if NULL
  if(cl != "numeric"){
    wtp <- mm.var$he$k
  }
  wtp.length <- length(wtp)

  # Format N
  cl <- class(N)
  # Select all N values if NULL
  if(cl != "numeric"){
    N.min <- min(mm.var$N.size)
    N.max <- max(mm.var$N.size)
    #if(N.min == N.max){
    #  N <- N.min
    #}
    N <- unique(round(seq(N.min, N.max, length.out=100)))
  }
  N.length<-length(N)
  
  N.length.calc<- length(unique(mm.var$N.size))
  
  # Format CI
    cl <- class(CI)
    # Set the default CI values if NULL
    if(cl != "numeric"){
      if(N.length.calc == 1){
        CI <- "No Uncertainty"
      }
      else{
        CI <- c(0.025, 0.25, 0.5, 0.75, 0.975)
      }
    }
    #Set N as requested by user
    if(cl == "numeric"){
      CI <- CI
      CI.inv <- 1 - CI
      CI <- ordered(unique(CI, CI.inv))
    }
    CI.length <- length(CI)

  ### Calcualte the EVSI
  e.var <- var(mm.var$he$delta.e)
  c.var <- var(mm.var$he$delta.c)
  e.fit <- as.matrix(mm.var$evppi$fitted.effects[, -mm.var$he$n.comparators])
  c.fit <- as.matrix(mm.var$evppi$fitted.costs[, -mm.var$he$n.comparators])
    
  simplify.var <- simplify2array(mm.var$variance.Q)
  
  if(N.length.calc == 1){
    ### One sample size
    mean.var <- apply(simplify.var, 1:2, mean)
    mean.var.e <- mean.var[1:mm.var$he$n.comparisons, 1:mm.var$he$n.comparisons]
    mean.var.c <- mean.var[(mm.var$he$n.comparisons + 1):(2 * mm.var$he$n.comparisons), (mm.var$he$n.comparisons + 1):(2 * mm.var$he$n.comparisons)]
    
    N.true <- N
    N <- NA
    index <- matrix(seq(1, mm.var$he$n.comparisons^2), nrow = mm.var$he$n.comparisons)
    }
  
  if(N.length.calc > 1){
  ## Multiple sample sizes
  n.unique.entries <- mm.var$he$n.comparisons * (mm.var$he$n.comparisons + 1) / 2
  beta.mat <- array(NA,dim = c(3000, 2, n.unique.entries))
  index <- matrix(NA, nrow = mm.var$he$n.comparisons, ncol = mm.var$he$n.comparisons)
  
  
  n.entry <- 1
  for(i in 1:mm.var$he$n.comparisons){
    for(j in i:mm.var$he$n.comparisons){
      beta.mat[,1,n.entry] <- fit.var.regression(
        simplify.var[1:mm.var$he$n.comparisons, 1:mm.var$he$n.comparisons,],
        e.var, e.fit, mm.var$N.size, i, j, mm.var$he, mm.var$update)
      beta.mat[,2,n.entry] <- fit.var.regression(
        simplify.var[(mm.var$he$n.comparisons + 1):(2 * mm.var$he$n.comparisons),
        (mm.var$he$n.comparisons + 1):(2 * mm.var$he$n.comparisons),], 
        c.var, c.fit, mm.var$N.size, i, j, mm.var$he, mm.var$update)
      # How to populate the variance matrices
      index[i, j] <- index[j, i] <- n.entry
      n.entry<-n.entry+1
    }
  }

  # Use CI to create a beta matrix
  # Find the appropriate quantiles for the beta parameter of interest
  beta.quantiles <- array(apply(beta.mat, c(2,3), quantile, prob=CI),
                          dim = c(CI.length, 2, n.unique.entries))

  #Extracting the beta parameter for the costs and effects
  mean.var.e<-beta.quantiles[,1,]
  mean.var.c<-beta.quantiles[,2,]
  }

  e.rescaled <- array(NA, dim=c(dim(e.fit), N.length, CI.length))
  c.rescaled <- array(NA, dim=c(dim(c.fit), N.length, CI.length))
  for(i in 1:CI.length){
    e.rescaled[, , , i] <- sapply(N,rescale.fitted,
                        mean.var = matrix(as.matrix(mean.var.e)[i, index],
                                          nrow = mm.var$he$n.comparisons,
                                          ncol = mm.var$he$n.comparisons),
                        fit = e.fit,
                        prior.var = e.var,
                        he = mm.var$he)
    c.rescaled[, , , i] <- sapply(N,rescale.fitted,
                                  mean.var = matrix(as.matrix(mean.var.c)[i, index],
                                                    nrow = mm.var$he$n.comparisons,
                                                    ncol = mm.var$he$n.comparisons),
                                  fit = c.fit,
                                  prior.var = c.var,
                                  he = mm.var$he)
  }


  wtp.func<-function(wtp.s,i,j,he){
    INB.star <- wtp.s * (- e.rescaled[, , j, i]) + c.rescaled[, , j, i]
    INB.augment <- as.data.frame(cbind(INB.star,0))
    mean.func <- function(x){ sum(x) / he$n.sim }
    EVSI <- sum(do.call(pmax, INB.augment))/he$n.sim -
      max(apply(INB.augment, 2, mean.func))
    return(EVSI)
  }

  EVSI.mat <- array(NA, dim = c(N.length, wtp.length, CI.length))
  for(i in 1:CI.length){
    for(j in 1:N.length){
        EVSI.mat[j, , i] <- sapply(wtp, wtp.func, i = i, j = j, he = mm.var$he)
    }
  }

  if(is.na(N[1])){
    N <- N.true
  }
  
  
  to.return<-list(evsi = EVSI.mat,
                  attrib = list(wtp = wtp, N = N, CI = CI),
                  evppi=mm.var$evppi,
                  he=mm.var$he)
  class(to.return)<-"evsi"
  return(to.return)
}

