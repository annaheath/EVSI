##summary.evsi################################################################
summary.evsi<-function(obj.evsi, wtp = NULL, N = NULL, CI = NULL){
  ##'Summarises the output of an evsi object
  ##INPUTS:
  ##'@param obj.evsi An evsi object.
  ##'@param wtp A fixed willingness-to-pay for which we would like to summarise the
  ##'  EVSI. If NULL then the willingness-to-pay assocaited with the optimal decision
  ##'  uncertainty is chosen.
  ##'@param N The sample size for which the EVSI is to be summerised. If NULL then chose
  ##'  the middle value for N.
  ##'
  ##OUTPUTS:
  ##' @return Prints a summary table with information on the Value of Information analysis
  ##' @author Anna Heath and Gianluca Baio
  ##' @export summary.evsi

  ## Check conditions
  if(is.null(wtp)){
    wtp <- obj.evsi$he$kstar[1]
  }
  
  if(!wtp %in% obj.evsi$attrib$wtp) {
    wtp.orig <- wtp
    wtp <- obj.evsi$attrib$wtp[which.min((obj.evsi$attrib$wtp-wtp)^2)]
    cat(paste("Note: The EVSI has not been computed for a willingness to pay of ", wtp.orig,
        " the VoI analysis will be displayed for a willingness to pay of ", wtp, "\n",sep=""))
  }
  
  if(is.null(N)){
    N <- quantile(obj.evsi$attrib$N, probs = 0.5, type = 4)
  }
  if(!N %in% obj.evsi$attrib$N) {
    N.orig <- N
    N <- obj.evsi$attrib$N[which.min((obj.evsi$attrib$N-N)^2)]
    cat(paste("Note: The EVSI has not been computed for a sample size of ", N.orig,
              " the EVSI analysis will be displayed for a sample size of ", N, "\n",sep=""))
  }
  
  if(!class(obj.evsi$attrib$CI) == "character"){
    if(is.null(CI)){
      CI <- obj.evsi$attrib$CI[which.min((obj.evsi$attrib$CI - 0.5)^2)]
    }
  }
  
  if(class(obj.evsi$attrib$CI) == "character"){
    CI <- obj.evsi$attrib$CI
  }
  
    if(!CI %in% obj.evsi$attrib$CI){
      CI <- obj.evsi$attrib$CI[which.min((obj.evsi$attrib$CI - 0.5)^2)]
    }

  ## Value of Information Selections
  evpi <- obj.evsi$he$evi[obj.evsi$he$k == wtp]
  evppi <- obj.evsi$evppi$evppi[obj.evsi$he$k == wtp]
  evsi <- obj.evsi$evsi[obj.evsi$attrib$N == N,obj.evsi$attrib$wtp == wtp, obj.evsi$attrib$CI == CI]
  
  ## Prints the summary table
  
  cat("\n")
  cat("Value of Information Summary\n")
  cat("\n")
  cat(paste("For willingness to pay parameter k = ",wtp,"\n",sep=""))
  cat("\n")
  cat(paste("For sample size of proposed trial N = ",N,"\n",sep=""))
  cat("\n")
  cat(paste("EVPI: ",evpi,"\n",sep=""))
  cat(paste("EVPPI: ",evppi,"\n",sep=""))
  cat(paste("EVSI: ",paste(evsi,"\n",sep="",collapse = ", "),sep = ""))
  
}
