data.quadrature <- function(index.data, data, sample, Q, N){
  ##'@param index.data - the column indices for the data in the prior preditive sample
  ##'@param data - character vector of the data names
  ##'@param sample - the dataframe from the Bayesian updating
  ##'@param Q - the number of posterior updates for the Heath et al method.
  ##'@param N - the sample size(s) of the future datasets
  
  # Number of data sets
  dat.num <- sapply(index.data, length)
  # Set up data objects
  data.full <- list()
  data.future <- vector("list",Q)
  
  if(length(index.data) > 1){
    matrix.data.test <- (length(grep(",",colnames(sample[,unlist(index.data)]))) != 0)
    cors <- sapply(1:min(dat.num),function(i){
      cor(sample[,index.data[[1]][i]],sample[,index.data[[2]][i]])
    })
    correlation.test <- sign(quantile(cors, probs = c(0.025))) == sign(quantile(cors, probs = c(0.975)))
    
    if(matrix.data.test | correlation.test){
      warning("Data generation within the EVSI package is not suitable for 
            correlated data, please use the gen.quantiles function to
            generate the data externally.")
    }
  }
  # Create a data matrix with only future data in it
  for(l in 1:length(data)){
  data.mat <- as.matrix(sample[, index.data[[l]]])
  categorical.check <- length(unique(as.numeric(as.matrix(data.mat))))
  
  if(categorical.check > 4){
    data.quants <- apply(data.mat, 2, quantile,probs = 1:Q/(Q+1),type = 4)
    cluster.centre <- array(NA, dim = c(Q, dat.num[l]))
    for(i in 1:dim(data.mat)[2]){
      cluster.centre[,i] <- sample(data.quants[,i],replace=FALSE)
    }

  }
  if(categorical.check <= 4){
    # Determine a cluster allocation for the data points
    cluster.allocation <- kmeans(data.mat, centers = Q, iter.max = 40)$cluster
    # Extract the centre of each cluster
    cluster.centre <- array(NA, dim = c(Q, sum(dat.num)))

    for(j in 1:Q){
      cluster.centre[j,] <- apply(data.mat[which(cluster.allocation == j), , drop = FALSE], 2, sample, 1)
    }
  }
    data.full[[l]] <- cluster.centre 
  }
  
  for(q in 1:Q){
    for(l in 1:length(data)){
      data.future[[q]][[l]]<-data.full[[l]][q,1:N[q]]
    }
    names(data.future[[q]])<-data
  }
  return(data.future)
}
