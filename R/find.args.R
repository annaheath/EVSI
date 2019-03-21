# Set the arguements of a function from the output of a Bayesian updating
find.args <- function(func, parameter.sims, monitor, row){
  ##'@param func - the function of interest
  ##'@param parameter.sims - a dataframe of parameter simulations from Bayesian updating
  ##'@param monitor.orig - a vector of the parameters that were monitored for Bayesain updating
  ##'@param row - the row number of the parameter simulations that contains the required function arguements
  
  # Function to extract matrix parameters from the parameter.sims dataframe
  extract.param.mat <- function(param.mat, param.names){
    # param.mat - a string with the name of a matrix parameter
    # param.names - a character vector with the names of all parameters
    param.values <- as.numeric(
      parameter.sims[row, grep(paste("\\b",param.mat,"\\b",sep=""), param.names)]
    )
    return(param.values)
  }
  
  # Find the function arguements
  args.names <- names(formals(func))
  param.names <- names(parameter.sims)
  
  # Determine which columns contain vector parameters, i.e. parameters that only appear once in the 
  # dataframe of parameter simulations
  param.single <- which(sapply(monitor, grep.fun, model.file.text = param.names) == 1)
  mat.single <- parameter.sims[row,monitor[param.single]]
  
  # Parameters that are matrices
  param.multi <- monitor[-param.single]
  
  #Sets the function arguements for the single parameters
  formals(func)[which(args.names %in% monitor[param.single])]<-
    mat.single[which(monitor[param.single] %in% args.names)]
  
  #Sets the function arguements for the matrix parameters
  multi.names <- args.names[args.names %in% param.multi]
  if(length(multi.names)!=0){
    for(i in 1:length(multi.names)){
      formals(func)[[multi.names[i]]] <-
        extract.param.mat(multi.names[i], param.names)
    }
  }
  return(formals(func))
}