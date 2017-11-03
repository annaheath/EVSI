##evsi.pop############################################################
evsi.pop<-function(evsi,trial.cost=NULL,setup=NULL,pp=NULL,
                   Pop=c(0,10000),Time=c(1,20),Dis=0.035,
                   wtp=NULL,N=NULL,pos=c("topright")){
  ##'Produces a graphic that considers the cost-effectiveness of a single future trial for different
  ##'Time Horizons for the treatment and incidence population
  ##'INPUTS:
  ##'@param evsi An evsi object.
  ##'@param trial.cost Either a fixed trial cost (numeric),a minimum and maximum value for the possible trial costs
  ##'    (numeric) or a distribution for a trial cost (character string) - to give probability of cost-effectiveness
  ##'    for the trial. Can be NULL if setup and pp are both given (and N if the sample size is not specified when
  ##'    the EVSI has been calculated/uploaded).
  ##'@param setup Gives the setup costs of the trial. Can be NULL if trial.costs are given. If NULL and pp given
  ##'   setup is taken to be 0.
  ##'@param pp Gives the per person costs of the trial. Can be NULL if trial.costs are given.
  ##'@param Pop The minimum and maximum possible values for the number of people benefitting from
  ##'   the treatment
  ##'@param Time The minimum and maximum possible values for the time horizon of the treatment.
  ##'@param Dis The discount rate for future studies, based on NICE guidelines.
  ##'@param wtp The willingness to pay value that this analysis is being undetaken for.
  ##'   If NULL then the function will automatically select the central wtp to default in BCEA will be 25000
  ##'@param N The sample size which we want to consider for the trial. If NULL then taken as the median trial size
  ##'   considered.
  ##'@param pos The position where the legend will be printed in the resulting graph.
  ##'
  ##OUTPUTS:
  ##' @return A graphic that gives the probability of a cost-effective trial.

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
  if(class(evsi$attrib$N)!="numeric"){
    N.select<-1
    if(class(N)!="numeric"){N<-evsi$attrib$N}

  }
  if(class(evsi$attrib$N)=="numeric"){
    if(class(N)!="numeric"){
      N.select<-ceiling(length(evsi$attrib$N)/2)
      N<-evsi$attrib$N[N.select]
    }
    if(class(N)=="numeric"){
      N.select<-which.min((evsi$attrib$N-N)^2)
    }
  }


  #Select deteministic or probabilistic EVSI
  if(length(evsi$attrib$CI)==1){
    type.evsi<-"det"
    ##Select per person EVSI for plot
    evsi.params<-c(evsi$evsi[N.select,wtp.select,1],0)}

  if(length(evsi$attrib$CI)>1){
    type.evsi<-"rand"
    q.1<-qnorm(evsi$attrib$CI[1])
    q.2<-qnorm(evsi$attrib$CI[2])
    evsi.1<-evsi$evsi[N.select,wtp.select,order(evsi$evsi[N.select,wtp.select,])[1]]
    evsi.2<-evsi$evsi[N.select,wtp.select,order(evsi$evsi[N.select,wtp.select,])[2]]
    evsi.params<-c((q.2*evsi.1-evsi.2*q.1)/(q.2-q.1),(evsi.2-evsi.1)/(q.2-q.1))
    #quants<-array()
    #for(i in 1:(length(evsi$attrib$CI)-1)){
    #  quants[i]<-pnorm(evsi$evsi[N.select,wtp.select,i],evsi.params[1],evsi.params[2])-evsi$attrib$CI[length(evsi$attrib$CI)-(i-1)]}
    #discrep<-mean(quants,na.rm=T)
    #if(discrep>0.2){warning("EVSI distribution not well approximated by a normal, the probabilistic EVSI analysis may be incorrect")}
  }

  ##Determine trial costs
  if(is.null(trial.cost)){
    if(class(N)=="character"){
      stop("Please define the trial costs using trial.costs or the sample size of experiment using N=")
    }

    if(is.null(setup)||is.null(pp)){stop("Please give the trial costs using either trial.costs for the full costs
                                         or setup and pp to give the set up and per person costs ")}

    setup.params<-c(mean(setup),(range(setup)[2]-range(setup)[1])/4)
    pp.params<-c(mean(pp),(range(pp)[2]-range(pp)[1])/4)
    trial.cost<-c(setup.params[1]+pp.params[1]*N,sqrt(setup.params[2]^2+N^2*pp.params[2]^2))
    }


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

  if(class(trial.cost)=="numeric"){
    if(length(trial.cost)==1){trial.cost.params<-c(trial.cost,0)}
    if(length(trial.cost)==2){trial.cost.params<-c(mean(trial.cost),(range(trial.cost)[2]-range(trial.cost)[1])/2)}

    Time.min<-min(Time);Time.max<-max(Time);Pop.min<-min(Pop);Pop.max<-max(Pop)
    dens.points<-100
    Time.seq<-seq(Time.min,Time.max,length.out=dens.points)
    Pop.seq<-seq(Pop.min,Pop.max,length.out=dens.points)

    Prob.mat<-matrix(NA,nrow=length(Time.seq),ncol=length(Pop.seq))
    for(i in 1:length(Time.seq)){
      for(j in 1:length(Pop.seq)){
        ENBS.mean<-ENBS(evsi.params[1],Pop.seq[j],Time.seq[i],Dis,trial.cost.params[1])
        ENBS.sd<-ENBS.sd.calc(evsi.params[2],Pop.seq[j],Time.seq[i],Dis,trial.cost.params[2])
        #Need to find the probability that the ENBS is GREATER than 0
        Prob.mat[i,j]<-pnorm(0,ENBS.mean,ENBS.sd,lower.tail=FALSE)
      }
    }

  }

  EVSI.full<-function(evsi,Pop,Time,Dis){
    evsi<-Pop*evsi/Dis*(1-exp(-Dis*Time))
    return(evsi)
  }
  if(class(trial.cost)=="character"){
    if(type.evsi=="rand"){
      warning("The trial cost given as a distribution so probablistic EVSI analysis cannot be performed")
      CI.select<-which.min((evsi$attrib$CI-0.5)^2)
      type.evsi<-"det"
      ##Select per person EVSI for plot
      evsi.focal<-evsi$evsi[N.select,wtp.select,CI.select]

    }

    split.1<-strsplit(trial.cost,split=NULL)
    commas<-which(split.1[[1]]==",")
    first.braket<-which(split.1[[1]]=="(")+1
    second.braket<-which(split.1[[1]]==")")-1
    n.args<-length(commas)+1
    args<-list()
    args[[1]]<-eval(parse(text=paste(split.1[[1]][first.braket:(commas[1]-1)],sep="",collapse="")))
    if(n.args>2){for(i in 2:(n.args-1)){
      args[[i]]<-(paste(split.1[[1]][(commas[i-1]+1):(commas[i]-1)],sep="",collapse=""))
    }
    }
    args[[n.args]]<-eval(parse(text=(paste(split.1[[1]][(commas[length(commas)]+1):second.braket],
                                           sep="",collapse=""))))

    func.samp<-paste(split.1[[1]][2:(first.braket-2)],sep="",collapse="")
    func.samp<-paste("p",func.samp,sep="",collapse="")

    Time.min<-min(Time);Time.max<-max(Time);Pop.min<-min(Pop);Pop.max<-max(Pop)
    dens.points<-100
    Time.seq<-seq(Time.min,Time.max,length.out=dens.points)
    Pop.seq<-seq(Pop.min,Pop.max,length.out=dens.points)

    Prob.mat<-matrix(NA,nrow=length(Time.seq),ncol=length(Pop.seq))
    for(i in 1:length(Time.seq)){
      for(j in 1:length(Pop.seq)){
        EVSI<-EVSI.full(evsi.focal,Pop.seq[j],Time.seq[i],Dis)
        Prob.mat[i,j]<-do.call(func.samp,append(as.list(EVSI),args))
      }
    }
  }
  #Plotting from Prob.mat
  image(x=Time.seq,y=Pop.seq,z=Prob.mat,col=colours,
        main="Probability of Cost-Effective Trial",xlab="Time Horizon",ylab="Incidence Population",
        xlim = c(Time.min,Time.max),ylim=c(Pop.min,Pop.max),
        breaks=seq(0,1,length.out=101))

  legend(alt.legend,c("Prob=0",rep(NA,98/2),"Prob=.5",rep(NA,96/2),"Prob=1"),fill=colours,border=colours,cex=0.75,y.intersp=0.15,
         box.lwd = 0,box.col = "white",bg = "white")
  box()

}
