gen.quantiles <- function(parameter, param.mat, Q, N.size = NULL){
  ##'Compute the posterior variance across Q difference proposed data collection exercises 
  ##'This can be used to estimate the EVSI using the Heath et al. method.
  ##INPUTS:
  ##'@param parameter A vector of parameters that define the sampling distribution of the future
  ##'  data. This can be given as a string (or vector of strings) of names or a numeric vector, 
  ##'  corresponding to the column numbers of parameters of interest.
  ##'@param param.mat A matrix containing the simulations for the parameters. The matrix should 
  ##'  have the column names matching the names of the parameters.
  ##'@param Q The number of quadrature points used the estimate the EVSI.
  ##'@param N.size If the EVSI is to be computed across sample size, N.size is a vector that 
  ##'  defines the maximum and minimum sample sizes for the proposed future study.
  ##'  If NULL then the the sample size of the future data will not be proposed by the gen.quantiles
  ##'  function.
  ##'
  ##'@return A dataframe containing the parameters simulations and sample size that can
  ##'  be used to generate future data to compute the EVSI using the Heath et al. moment
  ##'  matching EVSI computation method.
  ##'@examples NULL
  
  ### Check conditions
  if (is.null(colnames(param.mat))) {
    colnames(param.mat) <- paste0("theta",1:dim(param.mat)[2])
  }
  
  if (class(parameter[1]) == "numeric" | class(parameter[1]) == "integer") {
    parameters = colnames(param.mat)[parameter]
  }
  
  else {
    parameters = parameter
    for (i in 1:length(parameters)) {
      parameter[i] <- which(colnames(param.mat) == parameters[i])
    }
    class(parameter)<-"numeric"
  }

    quantiles.parameters <- array(NA, dim = c(Q, length(parameter)))
  colnames(quantiles.parameters) <- parameters
  for(i in 1:length(parameter)){
    quantiles.parameters[,i] <- sample(quantile(param.mat[,parameter[i]],
                                                probs = 1:Q / (Q + 1), type = 4))
  }
    
    if (!is.null(N.size)) {
      N.size <- trunc((seq(sqrt(min(N.size)),sqrt(max(N.size)),length.out = Q))^2)
      quantiles.aug <- cbind(quantiles.parameters, N = N.size)
      quantiles.parameters <- quantiles.aug
    }
  return(as.data.frame(quantiles.parameters))
    }