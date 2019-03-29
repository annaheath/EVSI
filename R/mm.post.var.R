##mm.post.var################################################################
mm.post.var <- function(model.stats, data, N.name = NULL, N.size = NULL, 
                        effects = NULL, costs = NULL, he = NULL, evi = NULL, parameters = NULL,
                        Q = 30, data.stats = NULL, update = c("bugs", "jags"),
                        n.burnin=1000, n.thin=1, n.iter=5000){
  ##'Compute the posterior variance across Q difference proposed data collection exercises 
  ##'This can be used to estimate the EVSI using the Heath et al. method.
  ##INPUTS:
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or JAGS model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. If the future data has already been generated then
  ##'   the data arguement can be given as a list. The elements of the data list should
  ##'   be data lists for JAGS or BUGS.
  ##'@param N.name A string that defines the name of the variable that defines the sample
  ##'    size of the future data in the model file. If NULL then the sample size of the data
  ##'    must be defined in data.stats.
  ##'@param N.size A scalar or vector that defines the sample sizes for the future study.
  ##'   A scalar N.size will compute the EVSI for a single sample size. If a vector is 
  ##'   passed for N.size then the EVSI can be estimated for all sample sizes within the
  ##'   range of values given. If NULL then the sample size of the data must be defined in
  ##'   data.stats.
  ##'@param effects This can either be given as a string which defines the name of the
  ##'   effectivness measure in the BUGS/JAGS model. Or it can be given as a
  ##'   function that takes the output of the model and calculates the
  ##'   effectiveness measure. The inputs of this function should be parameter
  ##'   names found in the model file.
  ##'@param costs This has a similar format to e but gives the value of the costs in
  ##'   the model file or a function defining the costs.
  ##'@param he A bcea object containing the base case analysis for the model.
  ##'@param evi An evppi object that contains the EVPPI analysis for the parameters
  ##'   that are being informed by the sample information.
  ##'@param parameters A list of the names of the parameters that are being targetted by the study.
  ##'  This can be given OR an evppi object needs to be given.
  ##'@param Q The number of quadrature points used the estimate the EVSI.
  ##'@param data.stats A data file for the BUGS/JAGS model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the JAGS/BUGS model
  ##'@param n.burnin The burnin for the JAGS/BUGS model
  ##'@param n.iter The number of interations for the JAGS/BUGS model
  ##'
  ##'@return An var.mm object.
  ##'1. evsi An array containing the EVSI by wtp, N and across different uncertaincies
  ##'2. attrib A list of wtp, N and prob describing the attributes of the evsi matrix.
  ##'3. evppi An evppi object containing all the information about the calculation of
  ##'   EVPPI.
  ##'4. he A bcea object containing all the information about the underlying health
  ##'   economic model
  ##'@examples NULL

  ### Check conditions
  
  # Future data provided by user.
  cl.dat<-class(data)
  if((cl.dat=="list")&&(length(data)!=Q)){
    warning(paste("The number of simulations for the future data does not equal the specified number of MC simulations
                  - the number of similations will be updated to ",length(data),".",sep=""))
    Q<-length(data)

  }

  # Functions for costs and effects not provided by user.
  if(class(costs)!="function"){
    costs<-function(costs){
      return(costs)
    }
  }

  if(class(effects)!="function"){
    effects<-function(effects){
      return(effects)
    }
  }
  
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
  ### Set up
  
  # Extract the arguements of the effects and costs functions
  argument.names <- unique(c(names(formals(effects)), names(formals(costs))))
  # Model parameters that are updated in Bayesian model
  monitor <- names(which(sapply(argument.names, grep.fun, model.file.text = readLines(model.stats)) > 0))
  
  # Set sample size if provided by user
  if(!is.null(N.name)){
    if(is.null(N.size)){"Please specify the required sample size estimation for the proposed data collection exercise."}
    N.size <- round(exp(seq(log(min(N.size)),log(max(N.size)),length.out = Q)))
    N.max <- list(max(N.size))
    names(N.max) <- N.name
    data.stats.update <- append(data.stats, N.max)
  }
  
  ### Prior predictive distribution

  # When prior predictive sampling is required
  if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
    if(is.null(N.size)){data.stats.update <- data.stats}
    monitor.data <- monitor
    
    # If data to be generated within function
    if(cl.dat=="character"){
      monitor.data <- unique(c(monitor.data, data))
    }
    
    # If EVPPI to be computed within function
    if(length(parameters)!=0){
      monitor.data <- unique(c(monitor.data,parameters))
    }
    
    sample.pp <- update.func(model.stats, data.stats.update, monitor.data, n.burnin, n.iter, n.thin)

  # When data is sampled from prior predictive distribution
  if(cl.dat=="character"){
      # Which columns contain the data
      index.data<-list()
      for(l in 1:length(data)){
        index.data[[l]] <- grep(data[l], colnames(sample.pp))
      }
      
      # Set sample size if not provided by user
      if(is.null(N.name)){
        N.size <- rep(length(index.data[[1]]),Q)
      }
      
      data <- data.quadrature(index.data, data, sample.pp, Q, N.size)
    }

    if(class(he)!="bcea"){
      # Checks if BCEA is installed (and if not, asks for it)
      if(!isTRUE(requireNamespace("BCEA",quietly = TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }
      # Calculate costs and effects for all simulations
      i <<- 1
      # Set function arguements
      formals(effects) <- find.args(effects, sample.pp, monitor, i)
      # Runs model with all arguements set as above
      model.effects <- effects()
      e <- c. <- matrix(NA, nrow=dim(sample.pp)[1], ncol=length(model.effects))
      for(i in 1:dim(sample.pp)[1]){
        i <<- i
        formals(effects) <- find.args(effects, sample.pp, monitor, i)
        formals(costs) <- find.args(costs, sample.pp, monitor, i)
        e[i,] <- eval(effects(), envir=parent.frame())
        c.[i,] <- eval(costs(), envir=parent.frame())
      }
      
      he <- BCEA::bcea(e, c.)
    }

    if(class(evi)!="evppi"){
      if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }

      # Find the columns that contain parameters of interest
      index.param <- list()
      for(l in 1:length(parameters)){
        index.param[[l]] <- grep(parameters[l], colnames(sample.pp))
      }
      param.columns <- unlist(index.param)
      evi <- BCEA::evppi(param.columns, sample.pp, he)
    }

}
    ### Calculate posterior variance by quadrature
    var.prepost<-list()
    start<-Sys.time()
    for(q in 1:Q){
      # Both data sets must be given to jags
      data.full<-append(data[[q]],data.stats)
      if(!is.null(N.name)){
        ss <- list(N.size[q])
        names(ss) <- N.name
        data.full <- append(data.full, ss)
      }
      # Sample from Bayesian updating
      sample.dat <- update.func(model.stats, data.full, monitor, n.burnin, n.iter, n.thin)
      # Use the JAGS output as inputs for the effects and costs functions
      incremental.benefit <- matrix(NA, nrow=dim(sample.dat)[1], ncol=2*he$n.comparisons)

      
      #Calculate costs and effects for all simulations
      for(i in 1:dim(sample.dat)[1]){
        i <<- i
        formals(effects) <- find.args(effects, sample.dat, monitor, i)
        formals(costs) <- find.args(costs, sample.dat, monitor, i)
        model.effects <- eval(effects(), envir=parent.frame())
        model.costs <- eval(costs(), envir=parent.frame())
        incremental.benefit[i,] <- cbind(model.effects[-he$ref] - model.effects[he$ref], model.costs[- he$ref] - model.costs[he$ref])
      }
      
      if(q == 1){
        end <- Sys.time()
      comp.time <- end - start
      print(paste(
        c("Model updating requires ", round(comp.time,2), " seconds. The EVSI will be calculated using ", Q, " model updates. The remaining computation time is around ",
          round(comp.time * (Q-1) / 60, 0), " minutes. The current time is ", strftime(Sys.time())),
                  sep="", collapse = ""))
      }
      print(paste("Update", q, "completed"))
      
      # Calculate the preposterior variance matrix for the costs and effects
      var.prepost[[q]] <- var(incremental.benefit)
    }
  
  #Return EVSI, plus evppi object and bcea object to plot EVSI plus attrib which fits in with later objects..
  list.return<-list(variance.Q = var.prepost,
                  N.size = N.size,
                  evppi = evi,
                  he = he,
                  update = update.func)
  class(list.return)<-"mm.var"
  return(list.return)
}

