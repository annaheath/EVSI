data.quadrature <- function(index.data, data, sample, Q, N){
  ##'@param index.data - the column indices for the data in the prior preditive sample
  ##'@param data - character vector of the data names
  ##'@param sample - the dataframe from the Bayesian updating
  ##'@param Q - the number of posterior updates for the Heath et al method.
  ##'@param N - the sample size(s) of the future datasets
  
  # Number of data sets
  dat.num <- sapply(index.data, length)
  
  # Create a data matrix with only future data in it
  data.mat <- sample[, unlist(index.data)]
  # Determine a cluster allocation for the data points
  cluster.allocation <- kmeans(data.mat, centers = Q, iter.max = 40)$cluster
  # Extract the centre of each cluster
  cluster.centre <- array(NA, dim = c(Q, sum(dat.num)))
  categorical.check <- length(unique(as.numeric(as.matrix(data.mat))))
  
  if(categorical.check > 3){
    for(l in 1:Q){
      cluster.centre[l,] <- apply(data.mat[which(cluster.allocation == l), ], 2, quantile, probs = 0.5, type = 1)
    }
  }
  if(categorical.check <= 3){
    for(l in 1:Q){
      cluster.centre[l,] <- apply(data.mat[which(cluster.allocation == l), ], 2, sample, 1)
    }
  }
  
  N.max <- max(N)
  # Set up data objects
  data.full <- list()
  data.future<-vector("list",Q)
  
  for(l in 1:length(data)){
    data.full[[l]] <- cluster.centre[, (((l - 1) * N.max) + 1):(l * N.max)] 
  }
  
  for(q in 1:Q){
    for(l in 1:length(data)){
      data.future[[q]][[l]]<-data.full[[l]][q,1:N[q]]
    }
    names(data.future[[q]])<-data
  }
  
  return(data.future)
}