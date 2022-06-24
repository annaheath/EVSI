rescale.fitted<-function(mean.var, prior.var, fit, he, N = NA){
  ##'@param mean.var - variance of posterior OR fitted beta value to estimate variance of posterior
  ##'@param if not already computed directly.
  ##'@param prior.var - variance of prior incremental net benefit
  ##'@param fit - fitted values from EVPPI calculation
  ##'@param he - BCEA object
  ##'@param N - sample size of future data.
  
  fit.var <- var(fit)
  
  if(!is.na(N)){
    #Variance for a specific N and beta
    mean.var <- prior.var - fit.var * (N / (N + mean.var))
  }
  
  # Rescale fitted values
  if(he$n_comparisons == 1){
    prepost.var <- max(0, prior.var - mean.var)
    rescaled <- (fit - mean(fit)) / sd(fit) * sqrt(prepost.var) + mean(fit)
  }
  
  if(he$n_comparisons>1){
    prepost.var <- prior.var - mean.var
    prepost.eigen.decomp <- base::eigen(prepost.var)
    # Matrix sqrt - from https://stat.ethz.ch/pipermail/r-help/2007-January/124147.html
    prepost.var.sqrt <- prepost.eigen.decomp$vectors %*% diag(sqrt(pmax(0, prepost.eigen.decomp$values))) %*% t(prepost.eigen.decomp$vectors)
    prior.mean <- matrix(rep(base::colMeans(fit), he$n.sim), nrow=he$n.sim, byrow=TRUE)
    
    # Matrix square root inverse
    prior.eigen.decomp <- eigen(fit.var)
    fit.var.sqrt.inv <- chol2inv(chol(prior.eigen.decomp$vectors %*% diag(sqrt(prior.eigen.decomp$values)) %*% t(prior.eigen.decomp$vectors)))
    
    # Rescale effects
    rescaled <-  (fit - prior.mean) %*% fit.var.sqrt.inv %*%
      prepost.var.sqrt + prior.mean
  }
  return(rescaled)
}
