# Standardised function to run Bayesian updating using jags and return dataframe
running.jags <- function(model, data, variable.names, n.burnin, n.iter, thin, inits = NULL){
  ##'@param model - file path for model
  ##'@param data - a list of the data to inform model
  ##'@param variable names - a vector of variable names to be tracked
  ##'@param n.burnin - the number of iterations to be used as burnin for the Markov chains
  ##'@param n.iter - the number of iterations for the Markov chains
  ##'@param thin - the thinning interval for the monitors
  ##'@param inits - initial values for the MCMC algoritm, not used by JAGS but 
  ##'  required for compatilibility
  
  # Create a jags model object
  model.create <- rjags::jags.model(model, data = data, quiet = TRUE)
  # Burnin chains
  update(model.create, n.burnin, progress.bar = "none")
  # Monitor results
  samples <- rjags::coda.samples(model.create, variable.names, n.iter = n.iter, n.thin = n.thin, progress.bar = "none")
  # Create dataframe for results
  samples.df <- as.data.frame(samples[[1]])
  
  return(samples.df)
}
