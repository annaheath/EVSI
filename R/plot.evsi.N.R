##plot.evsi.N###########################################################
plot.evsi.N<-function(evsi,wtp=NULL,pos=c("bottomright"),CI=NULL){
  ##'Calculating the EVSI for a specific WTP giving the uncertainty bands across the different
  ##'samples sizes
  ##INPUTS
  ##'@param evsi Output of the comp.evsi.N function
  ##'@param wtp The willingness to pay value that the graphic should be produced for - it will
  ##'   be chosen if wtp=NULL.
  ##'@param pos The position where the legend will be printed (default='bottomright')
  ##'@param CI The indexes that we would like to take from the CI in the evsi object.
  ##'
  ##OUTPUTS
  ##'@return EVSI The EVSI calculated for a specific wtp with uncertainty estimates.
  ##'@return A graphical representation of the uncertainty.

  alt.legend <- pos
  if (is.numeric(alt.legend) & length(alt.legend) == 2) {
    temp <- ""
    if (alt.legend[2] == 0)
      temp <- paste0(temp, "bottom")
    else if (alt.legend[2] != 0.5)
      temp <- paste0(temp, "top")
    if (alt.legend[1] == 1)
      temp <- paste0(temp, "right")
    else temp <- paste0(temp, "left")
    alt.legend <- temp
    if (length(grep("^((bottom|top)(left|right)|right)$",
                    temp)) == 0)
      alt.legend <- FALSE
  }
  if (is.logical(alt.legend)) {
    if (!alt.legend)
      alt.legend = "topright"
    else alt.legend = "topleft"
  }


  #Pick wtp threshold if not selected.
  if(class(wtp)!="numeric"){
    wtp.select<-which.min(abs(evsi$he$kstar-evsi$attrib$wtp))
    wtp<-evsi$attrib$wtp[which.min(abs(evsi$he$kstar-evsi$attrib$wtp))]
  }
  if(class(wtp)=="numeric"){
    wtp.select<-which.min(abs(wtp-evsi$attrib$wtp))
    wtp<-evsi$attrib$wtp[which.min(abs(wtp-evsi$attrib$wtp))]
  }

  if(class(CI)=="numeric"){
    CI.select<-CI
    CI<-evsi$attrib$CI[CI]

  }

  if(is.null(CI)){
    CI.select<-1:length(evsi$attrib$CI)
    CI<-evsi$attrib$CI
  }

  if(class(evsi$attrib$N)=="character"){
    stop("This plot gives the EVSI for increasing sample size. Do not use on a single design.")
  }

  CI.length<-length(CI)
  #Extracting the EVSI values for the wtp of interest
  EVSI<-array(NA,dim=c(length(evsi$attrib$N),1,CI.length))
  EVSI[,1,(1:CI.length)]<-rbind(evsi$evsi[,wtp.select,CI.select])

  #Set up the plot
  plot(1,1,ylim=c(min(EVSI)*0.95,max(EVSI)*1.05),xlim=c(min(evsi$attrib$N),max(evsi$attrib$N)),
       col="white",xlab=expression("Sample Size"),ylab="Per Person EVSI",oma=c(0,0,-1,0),main="Expected Value of Sample Information across Sample Size")

  if(CI.length%%2==1){
    lwd<-c(1:ceiling(CI.length/2),(ceiling(CI.length/2)-1):1,1)
    lty<-c(ceiling(CI.length/2):1,2:ceiling(CI.length/2),1)
  }
  if(CI.length%%2==0){
    lwd<-c(1:(CI.length/2),(CI.length/2):1,1)
    lty<-c((CI.length/2):1,2:(CI.length/2),1)
  }

  if(length(evsi$attrib$N)<15){
    for(l in 1:CI.length){
      points(evsi$attrib$N,EVSI[,,l],pch=19,
             lwd=lwd[l],lty=lty[l])
    }
  }

  if(length(evsi$attrib$N)>=15){
    for(l in 1:CI.length){
      points(evsi$attrib$N,EVSI[,,l],type="l",
             lwd=lwd[l],lty=lty[l])
    }
  }

  fitted.PI<--(wtp*evsi$evppi$fitted.e-evsi$evppi$fitted.c)
  abline(h=mean(apply(fitted.PI,1,max))-max(apply(fitted.PI,2,mean)),col="springgreen",lwd=lwd[CI.length+1],lty=lty[CI.length+1])

  legend(alt.legend,c(as.character(CI),"EVPPI"),
         col=c(rep("black",CI.length),"springgreen"),lwd=lwd,lty=lty,
         box.lwd = 0,box.col = "white",bg = "white")
  box()

}
