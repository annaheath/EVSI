##optim.ss#############################################################
optim.ss<-function(evsi,setup,pp,Pop,Time,wtp=NULL,Dis=0.035){
  ##'Calculates the optimal sample size for the future trial.
  ##INPUTS
  ##'@param evsi An evsi object that calculates the evsi by sample size
  ##'@param setup The setup costs of the trial can be given as a vector to represent
  ##'   uncertainty in the costs
  ##'@param pp The per person cost of the trial, can be given as a vector to represent
  ##'   uncertainty in these costs.
  ##'@param Pop The incidence population level
  ##'@param Time The time horizon.
  ##'@param wtp The willingness-to-pay value
  ##'@param Dis The Discount rate (default at 0.035)
  ##OUTPUTS
  ##'@return SS.max The optimal sample size.
  ##'@return ENBS The Expected Net Benefit of Sampling at the optimal sample size
  ##'@return SS.I The interval for which the ENBS is within 5% of the maximum ENBS
  ##'@return ENBS.I The ENBS values which correspond to 5% of the maximum ENBS - 
  ##'   used for plotting purposes.
  
  #evsi object?
  if(class(evsi)!="evsi"){stop("evsi must be in the evsi class - please create with either the evsi.calc or evsi.upload functions")}
  
  CI.select<-which.min((0.5-evsi$attrib$CI)^2)
  CI<-evsi$attrib$CI[CI.select]
  #Pick wtp threshold if not selected.
  if(class(wtp)!="numeric"){
    wtp.select<-which.min(abs(evsi$he$kstar-evsi$attrib$wtp))
    wtp<-evsi$attrib$wtp[which.min(abs(evsi$he$kstar-evsi$attrib$wtp))]
  }
  if(class(wtp)=="numeric"){
    wtp.select<-which.min(abs(wtp-evsi$attrib$wtp))
    wtp<-evsi$attrib$wtp[which.min(abs(wtp-evsi$attrib$wtp))]
  }
  
  EVSI<-evsi$evsi[,wtp.select,CI.select]
  if((length(setup)>1)||(length(pp)>1)){
    setup<-mean(setup)
    pp<-mean(pp)
  }
  
  ENBS<-Pop*EVSI/Dis*(1-exp(-Dis*Time))-setup-pp*evsi$attrib$N
  max.select<-which.max(ENBS)
  max.less<-max.select-1
  max.greater<-max.select+1
  if((max.less<1)|(max.greater>length(ENBS))){
    N.max<-evsi$attrib$N[max.select]
    ENBS.max<-ENBS[max.select]
    warning("Optimal sample size is at the limit of the considered values for N. An alternative sample size may be optimal,
            please consider alternative values of N in the evsi.calc function.")
  }
  else{
    N.fit<-evsi$attrib$N[c(max.less,max.select,max.greater)]
    N2<-N.fit^2
    ENBS.fit<-ENBS[c(max.less,max.select,max.greater)]
    model.optim<-lm(ENBS.fit~N.fit+N2)
    N.max<-round(-model.optim$coefficients[2]/(2*model.optim$coefficients[3]))
    ENBS.max<-predict(model.optim,list(N.fit=N.max,N2=N.max^2))}
  
  ##Creating a tolerance for the flat part of the graphic
  tol<-ENBS.max-abs(ENBS.max*0.05)
  limits<-which(ENBS>tol)
  N.range<-range(evsi$attrib$N[limits])
  ENBS.range<-ENBS[which(evsi$attrib$N %in% N.range)]
  
  
  return(list(SS.max=N.max,ENBS=ENBS.max,SS.I=N.range,ENBS.I=ENBS.range))
  
}


