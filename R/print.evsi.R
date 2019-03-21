##summary.evsi################################################################
print.evsi<-function(obj.evsi){
  ##'Summarises the output of an evsi object
  ##INPUTS:
  ##'@param obj.evsi An evsi object.
  ##'
  ##OUTPUTS:
  ##' @return Prints the Value of Information analysis
  ##' @author Anna Heath and Gianluca Baio
  ##' @export print.evsi
  
  ## Set conditions
  wtp <- obj.evsi$he$kstar[1]
  wtp <- obj.evsi$attrib$wtp[which.min((obj.evsi$attrib$wtp-wtp)^2)]
  N <- quantile(obj.evsi$attrib$N, probs = 0.5, type = 4)

    if(!class(obj.evsi$attrib$CI) == "character"){
      CI <- obj.evsi$attrib$CI[which.min((obj.evsi$attrib$CI - 0.5)^2)]
    }
  
  if(class(obj.evsi$attrib$CI) == "character"){
    CI <- obj.evsi$attrib$CI
  }

  ## Value of Information Selections
  evpi <- obj.evsi$he$evi[obj.evsi$he$k == wtp]
  evppi <- obj.evsi$evppi$evppi[obj.evsi$he$k == wtp]
  evsi <- obj.evsi$evsi[obj.evsi$attrib$N == N,obj.evsi$attrib$wtp == wtp, obj.evsi$attrib$CI == CI]
  
  ## Prints the summary table
  cat(paste("Willingness to pay of ",wtp,"\n",sep=""))
  cat("\n")
  cat(paste("Sample size of ",N,"\n",sep=""))
  cat("\n")
  cat(paste("EVPI: ",evpi,"\n",sep=""))
  cat(paste("EVPPI: ",evppi,"\n",sep=""))
  cat(paste("EVSI: ",paste(evsi,"\n",sep="",collapse = ", "),sep = ""))
  
}
