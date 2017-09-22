##enbs.plot############################################################
enbs.plot<-function(evsi,setup,pp,prob=NULL,
                    Pop=10000,Time=10,Dis=0.035,
                    wtp=NULL,N=NULL,pos=c("bottomright")){
  ##'Produces a graphic that considers the cost-effectiveness of a single future trial for different
  ##'Time Horizons for the treatment and incidence population
  ##'INPUTS:
  ##'@param evsi An evsi object.
  ##'@param setup Gives the minimum and maximum possible setup costs of the trial.
  ##'@param pp Gives the minimum and maxiumum possible per person costs of the trial.
  ##'@param prob The credible interval bounds that you would like plotted.
  ##'    By default the CI intervals.
  ##'@param Pop The number of people benefitting from the treatment.
  ##'@param Time The minimum and maximum possible values for the time horizon of the treatment.
  ##'@param Dis The discount rate for future studies, based on NICE guidelines.
  ##'@param wtp The willingness to pay value that this analysis is being undetaken for.
  ##'   If NULL then the function will automatically select the central wtp to default in BCEA will be 25000
  ##'@param N The sample sizes for which we want to consider the ENBS. If NULL taken as all the
  ##'   sample sizes that the EVSI has been calculated for.
  ##'@param pos The position where the legend will be printed in the resulting graph.
  ##'
  ##OUTPUTS:
  ##' @return A graphic that gives the ENBS for the range of N values given.

  #Legend
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


  #Select WTP
  if(class(wtp)!="numeric"){
    wtp.select<-ceiling(length(evsi$attrib$wtp)/2)
    wtp<-evsi$attrib$wtp[wtp.select]
  }
  if(class(wtp)=="numeric"){
    wtp.select<-which.min((evsi$attrib$wtp-wtp)^2)
  }

  #Select N
  if(class(N)!="numeric"){
    N<-evsi$attrib$N
  }

  if(class(N)=="character"){
    stop("Please define the sample size of your experiment using N=")
  }

  length.N<-length(N)

  N.select<-array(NA,dim=length.N)
  for(i in 1:length.N){
    N.select[i]<-which.min((evsi$attrib$N-N[i])^2)
  }

  #Select probabilities
  if(class(prob)!="numeric"){
    prob<-evsi$attrib$CI
  }

  length.prob<-length(prob)
  #Select deteministic or probabilistic EVSI
  if(length(evsi$attrib$CI)==1){
    type.evsi<-"det"
    ##Select per person EVSI for plot
    evsi.params<-cbind(evsi$evsi[N.select,wtp.select,1],0)}

  if(length(evsi$attrib$CI)>1){
    type.evsi<-"rand"
    q.1<-qnorm(evsi$attrib$CI[1])
    q.2<-qnorm(evsi$attrib$CI[2])
    evsi.params<-array(NA,dim=c(length.N,2))
    for(i in 1:length.N){
      evsi.1<-evsi$evsi[N.select[i],wtp.select,length(evsi$attrib$CI)]
      evsi.2<-evsi$evsi[N.select[i],wtp.select,length(evsi$attrib$CI)-1]
      evsi.params[i,]<-c((q.2*evsi.1-evsi.2*q.1)/(q.2-q.1),(evsi.2-evsi.1)/(q.2-q.1))}
  }


  ##Determine trial costs
  if(length(setup)!=2||length(pp)!=2){
    stop("Please give minimum and maximum values for the setup and per person trial costs.")}

  ##Setting up the stochastic values for the trial costs
  setup.params<-c(mean(setup),(range(setup)[2]-range(setup)[1])/4)
  pp.params<-c(mean(pp),(range(pp)[2]-range(pp)[1])/4)
  trial.cost<-array(NA,dim=c(length.N, 2))
  for(i in 1:length.N){
    trial.cost[i,]<-c(setup.params[1]+pp.params[1]*N[i],sqrt(setup.params[2]^2+N[i]^2*pp.params[2]^2))}



  #Colours for plotting
  colours<-colorRampPalette(c("black","navy","blue","skyblue","aliceblue","white"))(100)

  ENBS<-function(evsi,Pop,Time,Dis,cost){
    #This calculation comes from a discounting of (exp(-Dis*Time)) but then you need to integrate over all
    #time points - hence divide by the discount rate and then multiply by some weird factor.
    enbs<-evsi*Pop/Dis*(1-exp(-Dis*Time))-cost
    return(enbs)
  }
  ENBS.sd.calc<-function(evsi.sd,Pop,Time,Dis,cost.sd){
    var<- (Pop/Dis*(1-exp(-Dis*Time)))^2*evsi.sd^2+cost.sd^2
    return(sqrt(var))
  }

  ENBS.mat<-array(NA,dim=c(length.N,5))
  for(j in 1:length.N){
    ENBS.mean<-ENBS(evsi.params[j,1],Pop,Time,Dis,trial.cost[j,1])
    ENBS.sd<-ENBS.sd.calc(evsi.params[j,2],Pop,Time,Dis,trial.cost[j,2])
    ENBS.mat[j,]<-qnorm(prob,ENBS.mean,ENBS.sd)
  }
  #Plotting
  if(length.prob%%2==1){
    lwd<-c(1:ceiling(length.prob/2),(ceiling(length.prob/2)-1):1,1)
    lty<-c(ceiling(length.prob/2):1,2:ceiling(length.prob/2),1)
  }
  if(length.prob%%2==0){
    lwd<-c(1:(length.prob/2),(length.prob/2):1,1)
    lty<-c((length.prob/2):1,2:(length.prob/2),1)
  }
  plot.new()

  plot.window(xlim=c(min(N),max(N)),ylim=c(min(ENBS.mat),max(ENBS.mat)))
  title(main="Expected Net Benefit of Sampling by Sample Size",xlab="Sample Size",ylab="ENBS")
  axis(side=2)

  if(length.N<15){
    for(l in 1:length.prob){
      points(N,ENBS.mat[,l],pch=19,
             lwd=lwd[l],lty=lty[l])
    }
  }

  if(length.N>=15){
    for(l in 1:length.prob){
      points(N,ENBS.mat[,l],type="l",
             lwd=lwd[l],lty=lty[l])
    }
  }

  abline(h=0,col="springgreen",lwd=lwd[length.prob+1],lty=lty[length.prob+1])

  legend(alt.legend,c(as.character(prob),"ENBS=0"),
         col=c(rep("black",length.prob),"springgreen"),lwd=lwd,lty=lty,
         box.lwd = 0,box.col = "white",bg = "white")
  box()
  optimal<-optim.ss(evsi,setup,pp,Pop,Time,wtp=wtp,Dis=Dis)
  axis(side=1)#,at=c(optimal$SS.I,optimal$SS.max,min(N),max(N)))

  a<-0.04
  poi<-(1+a)*min(ENBS.mat)-a*max(ENBS.mat)
  points(c(optimal$SS.I),c(poi,poi),type="l",col="red",lwd=3)
  points(c(optimal$SS.max),min(poi),pch=9,col="red",lwd=3)
  #points(c(optimal$SS.I),c(min(ENBS.mat),min(ENBS.mat)),type="l",col="red",lwd=3)
  #points(c(optimal$SS.max),min(ENBS.mat),pch=9,col="red",lwd=3)

}

