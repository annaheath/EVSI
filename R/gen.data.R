gen.data<-function(model.stats, data, N.name = NULL, N.size = NULL, moment.matching = TRUE,
                   Q = 30, data.stats = NULL, update = c("bugs", "jags"),
                   n.burnin=1000, n.thin=1, n.iter=5000){
  ##'Compute the posterior variance across Q difference proposed data collection exercises 
  ##'This can be used to estimate the EVSI using the Heath et al. method.
  ##INPUTS:
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or JAGS model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. 
  ##'@param N.name A string that defines the name of the variable that defines the sample
  ##'    size of the future data in the model file. If NULL then the sample size of the data
  ##'    must be defined in data.stats.
  ##'@param N.size A scalar or vector that defines the sample sizes for the future study.
  ##'   A scalar N.size will compute the EVSI for a single sample size. If a vector is 
  ##'   passed for N.size then the EVSI can be estimated for all sample sizes within the
  ##'   range of values given. If NULL then the sample size of the data must be defined in
  ##'   data.stats.
  ##'@param moment.matching If TRUE then the returned data can be used as inputs to the Heath
  ##'  et al. moment matching method for EVSI computation. If FALSE then a data.frame of the
  ##'  prior predictive distribution of the data will be returned.
  ##'@param Q The number of quadrature points used the estimate the EVSI.
  ##'@param data.stats A data file for the BUGS/JAGS model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the JAGS/BUGS model
  ##'@param n.burnin The burnin for the JAGS/BUGS model
  ##'@param n.iter The number of interations for the JAGS/BUGS model
  ##'
  ##'@return A list of future data object.
  ##'1. data A data frame or list of the simulated data for the proposed study.
  ##'2. N.calc A vector containing the sample sizes for the which the EVSI has 
  ##'   been computed.
  ##'@examples NULL
  
  ### Check conditions
  
  # Future data provided by user.
  cl.dat<-class(data)
  if(cl.dat!="character"){stop("Please define the names of the data variables in your model.")}
  
  # Bayesian updating software available
  if(update=="jags"){ 
    # Checks if rjags is installed (and if not, asks for it)
    if(!isTRUE(requireNamespace("rjags",quietly=TRUE))) {
      stop("You need to install the R package 'rjags' and the software 'JAGS'. \nPlease see http://mcmc-jags.sourceforge.net/ for instructions
           on installing 'JAGS' and then run in your R terminal:\n install.packages('rjags')")
    }
  }
  
  if(update=="bugs"){
    # Checks if R2OpenBUGS is installed (and if not, asks for it)
    if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
      stop("You need to install the R package 'R2OpenBUGS' and the software 'OpenBUGS'. \nPlease see http://www.openbugs.net/w/FrontPage for instructions
           on installing 'OpenBUGS' and then run in your R terminal:\n install.packages('R2OpenBUGS')")
    }
  }
  
  update.func <- match.fun(paste("running.", update, sep=""))
  
  # Set sample size if provided by user
  if(!is.null(N.name)){
    if(is.null(N.size)){stop("Please specify the required sample size estimation for the proposed data collection exercise.")}
    N.size <- trunc((seq(sqrt(min(N.size)),sqrt(max(N.size)),length.out = Q))^2)
    N.max <- list(max(N.size))
    names(N.max) <- N.name
    data.stats.update <- append(data.stats, N.max)
  }
  
  if(is.null(N.size)){data.stats.update <- data.stats}
  sample.pp <- update.func(model.stats, data.stats.update, data, n.burnin, n.iter, n.thin)
  
  # Which columns contain the data
  index.data<-list()
  for(l in 1:length(data)){
    index.data[[l]] <- grep(data[l], colnames(sample.pp))
  }
  
  # Set sample size if not provided by user
  if(is.null(N.name)){
    N.size <- rep(length(index.data[[1]]),Q)
  }
  
  if(moment.matching){data <- data.quadrature(index.data, data, sample.pp, Q, N.size)}
  if(!moment.matching){data <- sample.pp[unlist(index.data)]}
  
  return(list(data = data,
              N = N.size))
}