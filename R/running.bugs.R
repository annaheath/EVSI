# Standardised function to run Bayesian updating using BUGS and return dataframe
running.bugs<-function(model, data, variable.names, n.burnin, n.iter, thin){
  ##'@param model - file path for model
  ##'@param data - a list of the data to inform model
  ##'@param variable names - a vector of variable names to be tracked
  ##'@param n.burnin - the number of iterations to be used as burnin for the Markov chains
  ##'@param n.iter - the number of iterations for the Markov chains
  ##'@param thin - the thinning interval for the monitors
  
  # Monitor results from bugs updating
  samples <- R2OpenBUGS::bugs(data, 
                              inits = NULL,
                              parameters.to.save=variable.names,
                              model.file=model, 
                              n.burnin = n.burnin,
                              n.iter = n.iter+n.burnin, 
                              n.thin = thin,
                              n.chain = 1,
                              DIC = FALSE, 
                              debug = FALSE)
  # Create dataframe for results
  samples.df <- as.data.frame(samples$sims.array[,1,])
  
  return(samples.df)
}
