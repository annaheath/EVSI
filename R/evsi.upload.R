##evsi.upload#########################################################
evsi.upload<-function(e,c,parameter,input,EVSI.mat=NULL,model.stats=NULL,wtp=NULL,N=NULL){
  ##'Uploads EVSI either from a matrix in R or a .csv file.
  ##'Allows EVSI calculated outside of package to be uploaded so graphics can be used.
  ##'INPUTS
  ##'@param e A matrix containing simuations for the variable of clinical effectiveness for
  ##'   each inteventions being considered.
  ##'@param c A matrix containing simuations for the variable of cost for each inteventions
  ##'   being considered.
  ##'@param parameter A vector of parameters that are updated by the future dataset. This can
  ##'   be given as a string of names or a numeric vector.
  ##'@param input A matrix of PSA simulations from the underlying health economic model.
  ##'@param EVSI.mat A matrix of EVSI values. The COLUMNS relate the different willingness
  ##'   to pay values and the ROWS relate the different sample sizes.
  ##'@param model.stats A .csv file containing the matrix of EVSI values. Again the COLUMNS
  ##'   should relate the different WTP values and the ROWS to different sample sizes.
  ##'@param wtp A vector of willingness to pay values for which the EVSI was calculated.
  ##'   If NULL then the first row of the EVSI matrix is assumed to contain the WTP values.
  ##'@param N A vector of sample sizes for which the EVSI is calculated. If NULL then the
  ##'   first column of the EVSI matrix is assumed to contain the sample sizes.
  ##'
  ##OUTPUTS
  ##' @return An evsi object.
  ##'1. evsi An array containing the EVSI by wtp, N and across different uncertaincies
  ##'2. attrib A list of wtp, N and prob describing the attributes of the evsi matrix.
  ##'3. evppi An evppi object containing all the information about the calculation of
  ##'   EVPPI.
  ##'4. he A bcea object containing all the information about the underlying health
  ##'   economic model
  ##'   @examples
  ##'   ...

  if(is.null(EVSI.mat)&&is.null(model.stats)){
    stop("You have not provided the EVSI. Please either use EVSI.mat or model.stats to give a matrix of EVSI values")}

  #Read in EVSI matrix from file
  if(is.null(EVSI.mat)){
    EVSI.mat<-as.matrix(read.csv(model.stats,header=FALSE))
  }

  #Set wtp values if not given
  if(is.null(wtp)){
    wtp<-as.numeric(EVSI.mat[1,])
    EVSI.mat<-EVSI.mat[-1,]
    if(is.null(N)){
      wtp<-wtp[-1]
    }
  }

  #Set N values if not given
  if(is.null(N)){
    N<-EVSI.mat[,1]
    EVSI.mat<-EVSI.mat[,-1]
  }

  #Set EVSI as an array
  EVSI.arr<-array(EVSI.mat,c(dim(EVSI.mat),1))

  #Find BCEA object
  if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
    stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
  }
  he<-BCEA::bcea(e,c)

  #Calculate EVPPI object
  evi<-BCEA::evppi(parameter,input,he)

  to.return<-list(evsi=EVSI.arr,attrib=list(wtp=wtp,N=N,CI=0.5),evppi=evi,he=he)
  class(to.return)<-"evsi"
  return(to.return)
}
