##plot.evsi###########################################################
plot.evsi<-function(evsi,pos=c(0,0.8),N=NULL){
  ##'Plots the EVSI across different willingness-to-pay thresholds for a single design.
  ##INPUTS
  ##'@param evsi An evsi object that calculates the EVSI for a single design using the evsi function
  ##'@param pos An indcation where the legend should be printed
  ##'@param N A vector of values for which the EVSI should be plotted.
  ##OTHER GRAPHICAL INPUTS TO BE ADDED.
  ##'
  ##OUTPUTS
  ##'@return Produces a graphic containing the EVPI, the EVPPI for the key parameters and the EVSI for the
  ##'fixed design.
  
  if(class(evsi)!="evsi"){stop("plot.evsi must be used with an evsi object.")}
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
  
  N.length<-length(evsi$attrib$N)
  
  #Which EVSIs to plot
  if(class(N)=="numeric"){
    select.length<-length(N)
    select<-rep(NA,select.length)
    for(i in 1:select.length){
      select[i]<-which.min((evsi$attrib$N-N[i])^2)
    }
    warning("The EVSI is only calculated for sample sizes listed in evsi$attrib$N. If N is not in this list then EVSI calculated for closest possible sample size.")
    N<-evsi$attrib$N[select]
  }
  
  #Determine which EVSI curves to pick if N is NULL (all)
  if(is.null(N)){
    select<-1:N.length
    N<-evsi$attrib$N
  }
  
  select.length<-length(select)
  
  if(class(evsi$attrib$CI)=="character"){CI.select<-1}
  if(class(evsi$attrib$CI)=="numeric"){CI.select<-which.min((evsi$attrib$CI-0.5)^2)}
  #Select EVSI values to plot
  EVSI<-rbind(evsi$evsi[select,,CI.select])
  
  #Plot the EVPI and create graphic
  plot(evsi$he$k, evsi$he$evi, t = "l", xlab = "Willingness to pay",
       ylab = "", main = "Expected Value of Sample Information",
       lwd = 2, ylim = range(range(evsi$he$evi), range(evsi$evppi$evppi),range(EVSI)))
  
  #Plot the EVPPI on the graphic
  col = "black"
  points(evsi$evppi$k, evsi$evppi$evppi, t = "l", col = col, lty = 1)
  
  #Choose colours for the different sample sizes
  colours<-colorRampPalette(c("skyblue","blue","darkblue"))(select.length)
  
  if(length(evsi$attrib$wtp)<30){
    for(s in 1:select.length){
      points(evsi$attrib$wtp,EVSI[s,],pch=19,col=colours[s])
    }
  }
  if(length(evsi$attrib$wtp)>=30){
    for(s in 1:select.length){
      points(evsi$attrib$wtp,EVSI[s,],type="l",col=colours[s])
    }
  }
  
  if(select.length==1){  legend(alt.legend, c("EVPI", "EVPPI for focal parameters",paste("EVSI for sample size of",N)),
                                col = c("black", "black",colours),
                                cex = 0.7,box.lwd = 0,box.col = "white",bg = "white", lty = c(1, 1,1),
                                lwd = c(2, 1))}
  if(select.length>1){legend(alt.legend,legend=c(min(N),rep(NA,max(0,select.length-2)),max(N)),
                             fill=colours,border=colours,cex=0.75,
                             y.intersp=max(0.1,1.2/select.length),
                             box.lwd = 0,box.col = "white",bg = "white")}
  box()
}
