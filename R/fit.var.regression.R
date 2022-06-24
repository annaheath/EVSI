fit.var.regression <- function(mm.var, var.prior, fit, N.calc, row, column, he, update.func){
  ##'@param mm.var - an array containing the posterior variance 
  ##'@param var.prior - a scalar or matrix containing the variance of the prior
  ##'@param fit.var - the fitted values from evppi function
  ##'@param N.calc - the sample sizes for which the future data was estiamted
  ##'@param row - the row of the variance matrix for which we are estimating the variance/covariance
  ##'@param column - the column of the variance matrix for which we are estiamting the variance/covarience
  ##'@param he - the health economic model output from bcea
  ##'@param update.func - the function that will be used to do the Bayesian updating
  
  # Variance of fitted values
  var.fit <- var(fit)
  
  # Specify data for non-linear regression fit.
  var.PI=x=NULL
  if(he$n_comparisons>1){
    pre.var <- mm.var
    prepost.var <- t(var.prior[row, column] - pre.var[row, column, ])
  }
  if(he$n_comparisons==1){
    pre.var <- mm.var
    prepost.var <- rep(var.prior, length(pre.var)) - pre.var
  }
  
  if(sd(prepost.var)==0){
    beta.ab <- array(1E10,dim = c(3000,1))
  }
  
  if(sd(prepost.var)!=0){
    sigma.tau=1/sd(prepost.var)
  data.a.b <- list(sigma.mu = sd(prepost.var)/2,
                   sigma.tau = sigma.tau,
                   N = length(N.calc),
                   shape.Nmax = 0.05 / max(N.calc),
                   var.PI = as.numeric(as.matrix(var.fit)[row, column]),
                   Nmax.par = max(N.calc)/2,
                   y = as.vector(prepost.var),
                   x = as.vector(N.calc)
  )
  
  ## Specify Bayesian model
  model.ab <- c("model
                {
                beta ~ dnorm(Nmax.par, shape.Nmax)  I(0.00000E+00, )
                for (i in 1:N) {
                y[i] ~ dnorm(mu[i], tau)
                mu[i] <- var.PI * (x[i]/(x[i] + beta))
                }
                sigma ~ dt(sigma.mu, sigma.tau, 3)  I(0.00000E+00, )
                tau <- 1/(sigma*sigma)
                }
                ")
  file.curve.fitting <- file.path("model_curve_fitting_EVSI.txt")
  writeLines(model.ab,file.curve.fitting)
  
  ## Model set up
  n.burnin <- 1000  # Number of burn in iterations
  n.thin <- 1
  n.iter <- 3000 # Number of iterations per chain
  
  beta.ab <- update.func(file.curve.fitting, data.a.b, c("beta"), n.burnin,
                         n.iter, n.thin , 
                         inits = list(c(beta = 1, sigma = sd(prepost.var)/2))) 
  }
  
  return(beta.ab[,1])
}
