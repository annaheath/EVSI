#####################################################################
#Suite of functions to calculate the EVSI using Heath et al. methods#
#####################################################################

##v1.0 June 2016 AH
## (C) Anna Heath + contributions from ...

##Functions Included#################################################

##evsi
##plot.evsi
##evsi.pop
##comp.evsi.N
##evsi.calc
##plot.evsi.N

##evsi################################################################
library(foreach)
library(BCEA)
library(rjags)
library(R2OpenBUGS)
library(shiny)
library(shinythemes)

evsi<-function(model.stats,data,effects=NULL,costs=NULL,he=NULL,evi=NULL,parameters=NULL,Q=30,data.stats=NULL,update=c("BUGS","jags"),
               n.burnin=1000,n.thin=1,n.iter=5000){
  ##'Compute the EVSI for a fixed design using the Heath et al. method
  ##'INPUTS:
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or jags model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. If the future data has already been generated then
  ##'   the data arguement can be given as a list. The elements of the data list should
  ##'   be data lists for jags or BUGS.
  ##'@param effects This can either be given as a string which defines the name of the
  ##'   effectivness measure in the BUGS/jags model. Or it can be given as a
  ##'   function that takes the output of the model and calculates the
  ##'   effectiveness measure. The inputs of this function should be parameter
  ##'   names found in the model file.
  ##'@param costs This has a similar format to e but gives the value of the costs in
  ##'   the model file or a function defining the costs.
  ##'@param he A bcea object containing the base case analysis for the model.
  ##'@param evi An evppi object that contains the EVPPI analysis for the parameters
  ##'   that are being informed by the sample information.
  ##'@param parameters A list of the names of the parameters that are being targetted by the study.
  ##'  This can be given OR an evppi object needs to be given.
  ##'@param Q The number of quadrature points used the estimate the EVSI.
  ##'@param data.stats A data file for the BUGS/jags model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the jags/BUGS model
  ##'@param n.burnin The burnin for the jags/BUGS model
  ##'@param n.iter The number of interations for the jags/BUGS model
  ##'
  ##OUTPUTS:
  ##'@return The EVSI for the specific design by willingness-to-pay. This can be plotted.
  ##'@examples
  ##'...

  #Is data in the correct format??
  cl.dat<-class(data)
  if((cl.dat=="list")&&(length(data)!=Q)){
    warning(paste("The number of simulations for the future data does not equal the specified number of MC simulations
                  - the number of similations will be updated to",length(data),"."))
    Q<-length(data)

  }

  #Set the functions for costs and effects
  if(class(costs)!="function"){
    costs<-function(costs){
      return(costs)
    }
  }

  if(class(effects)!="function"){
    effects<-function(effects){
      return(effects)
    }
  }

  #Extract the arguements of the effects and costs functions
  nameofinterest<-unique(c(names(formals(effects)),names(formals(costs))))

  #Check which of these parameters are also in the model file
  #grep reads the model file finds the occurance of a string.
  grep.fun<-function(name,main){
    length(grep(name, main, value = FALSE))
  }

  #Set the function arguements for the model functions
  find.args<-function(func,sample,i){
    #Find the names of all the function arguements
    args.names<-names(formals(func))
    #Find the number of times a variable is containted in the jags model
    #This allows us to determine matrix contructions for certain variables.
    #The variable a vector - i.e. it only appears once in the jags output
    jags.params<-which(sapply(moniter,grep.fun,main=names(sample))==1)
    inputs<-sample[,moniter[jags.params]]

    #Matrix Parameters
    multiple.params<-moniter[-jags.params]

    #Function to extract the matrix parameters from the jags object
    multiple.params.extract<-function(multi.params){
      if(length(multi.params)==0){return(NULL)}
      else{
        list.params<-list()
        length(list.params)<-length(multi.params)
        for(j in 1:length(multiple.params)){
          list.params[j]<-list(as.numeric(sample[i,grep(multi.params[j],names(sample))]))
        }
        return(list.params)}
    }

    #Sets the effects function arguements for the single parameters
    formals(func)[which(args.names %in% names(inputs))]<-
      inputs[i,which(names(inputs) %in% args.names)]
    #Sets the effects function arguements for the matrix parameters
    formals(func)[args.names %in% multiple.params]<-
      multiple.params.extract(args.names[which(args.names %in% multiple.params)])
    return(formals(func))
  }

  #Generate future samples by finding prior-predictive distribution
  if(update=="jags"){#Model can be written in jags.
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
      #Runs the jags model based on the rjags package
      Model.JAGS<- jags.model(model.stats,data=data.stats,quiet=TRUE)
      update(Model.JAGS,n.burnin,progress.bar="none")
      Prior.Pred <- coda.samples(Model.JAGS, prior.pred.data, n.iter=n.iter,n.thin=n.thin,progress.bar="none")
      PP.sample<-as.data.frame(Prior.Pred[[1]])}

    if(cl.dat=="character"){
      #Determine which columns contain the data
      index.data<-list()
      for(l in 1:length(data)){
        index.data[[l]]<-grep(data[l],colnames(PP.sample))
      }

      #Finds the quadrature points and restructures the data
      Ord<-order(PP.sample[,unlist(index.data)[1]])
      Quad<-trunc(sample((1:Q)/(Q+1),replace=F)*length(Ord))
      #ordered contains the future data samples that will be used for the quadrature
      ordered<-as.matrix(PP.sample[Ord[Quad],])
    }

    if(class(he)!="bcea"){
      #Calculate costs and effects for all simulations
      i<-1
      #Set function arguements
      formals(effects)<-find.args(effects,PP.sample,i)
      #Runs model with all arguements set as above
      model.effects<-effects()
      e<-c<-matrix(NA,nrow=dim(PP.sample)[1],ncol=length(model.effects))
      for(i in 1:dim(PP.sample)[1]){
        #Set function arguements
        formals(effects)<-find.args(effects,PP.sample,i)
        #Runs model with all arguements set as above
        e[i,]<-effects()
        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,PP.sample,i)
        c[i,]<-costs()
      }
      he<-bcea(e,c)
    }
    #Check number of interventions
    if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
                                multi-decision problems")}

    if(class(evi)!="evppi"){
      #Determine which model rows contain statements defining the data
      #lines.imp<-list()
      #for(l in 1:length(data)){
      #  lines.imp[[l]]<-grep(data[l],readLines(model.stats))
      #}

      #Determine which parameters are related to the data in those rows.
      #lines.import<-readLines(model.stats)[unique(unlist(lines.imp))]
      #params.pi<-names(which(sapply(unique(nameofinterest),grep.fun,main=lines.import)>0))

      #Find the columns that contain these parameters
      index<-list()
      for(l in 1:length(parameters)){
        index[[l]]<-grep(parameters[l],colnames(PP.sample))
      }
      of.interest<-unlist(index)
      #PSA.mat<-as.matrix(PP.sample[,unlist(index)])
      evi<-evppi(of.interest,PP.sample,he)
    }


    ####Calculate EVSI by Quadrature####
    var.prepost<-list()
    start<-Sys.time()

    #foreach(q=1:Q,.packages="rjags",.export=ls(.GlobalEnv)) %do% {
    for(q in 1:Q){

      if(exists("ordered",mode = "numeric")){#Creating list of the future data to give to jags
        length.data<-dim(ordered)[2]
        names.data<-data#colnames(ordered)

        Data.Fut<-list()
        for(d in 1:length(index.data)){
          Data.Fut[[d]]<-as.numeric(ordered[q,index.data[[d]]])
        }
        names(Data.Fut)<-names.data
      }
      if(cl.dat=="list"){
        Data.Fut<-data[[q]]
      }

      #Both data sets must be given to jags
      data.full<-append(Data.Fut,data.stats)

      Model.JAGS<- jags.model(model.stats,data =  data.full,quiet = TRUE)
      update(Model.JAGS,n.burnin,progress.bar="none")
      samples <- coda.samples(Model.JAGS, moniter, n.iter=n.iter,n.thin=n.thin,n.chain=1,progress.bar="none")
      #Create a dataframe containing all the JAGS samples
      sample<-as.data.frame(samples[[1]])

      #Use the JAGS output as inputs for the effects and costs functions
      int<-matrix(NA,nrow=dim(sample)[1],ncol=2*he$n.comparisons)

      #Calculate costs and effects for all simulations
      #int<-foreach(i=1:dim(sample)[1],.combine="rbind") %do% {
      for(i in 1:dim(sample)[1]){
       #Set function arguements
        formals(effects)<-find.args(effects,sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        #Calculate the incremental effects.
        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,sample,i)
        model.costs<-costs()
        int[i,]<-cbind(model.effects[-he$ref]-model.effects[he$ref],model.costs[-he$ref]-model.costs[he$ref])
      }
      #Calculate the preposterior variance matrix for the costs and effects to calculate EVSI by WTP
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))   }
      print(paste("Update",q,"completed"))
      var.prepost[[q]]<-var(int)
    }

    }

  #Analysis in BUGS
  if(update=="BUGS"){#Model can be written in BUGS.
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
      #Runs the BUGS model based on the R2OpenBUGS package
      Model.BUGS<- bugs(data.stats,inits=NULL,parameters.to.save=prior.pred.data,
                        model.file=model.stats, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chain=1,
                        DIC=FALSE,debug=FALSE)
      PP.sample<-as.data.frame(Model.BUGS$sims.matrix)}


    if(cl.dat=="character"){
      #Determine which columns contain the data
      index.data<-list()
      for(l in 1:length(data)){
        index.data[[l]]<-grep(data[l],colnames(PP.sample))
      }

      #Finds the quadrature points and restructures the data
      Ord<-order(PP.sample[,unlist(index.data)[1]])
      Quad<-trunc(sample((1:Q)/(Q+1),replace=F)*length(Ord))
      #ordered contains the future data samples that will be used for the quadrature
      ordered<-as.matrix(PP.sample[Ord[Quad],])
    }

    if(class(he)!="bcea"){
      #Calculate costs and effects for all simulations
      i<-1
      #Set function arguements
      formals(effects)<-find.args(effects,PP.sample,i)
      #Runs model with all arguements set as above
      model.effects<-effects()
      e<-c<-matrix(NA,nrow=dim(PP.sample)[1],ncol=length(model.effects))
      for(i in 1:dim(PP.sample)[1]){
        #Set function arguements
        formals(effects)<-find.args(effects,PP.sample,i)
        #Runs model with all arguements set as above
        e[i,]<-effects()
        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,PP.sample,i)
        c[i,]<-costs()
      }
      he<-bcea(e,c)
    }
    #Check number of interventions
    if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
                                multi-decision problems")}

    if(class(evi)!="evppi"){
      #Determine which model rows contain statements defining the data
      #lines.imp<-list()
      #for(l in 1:length(data)){
      #  lines.imp[[l]]<-grep(data[l],readLines(model.stats))
      #}

      #Determine which parameters are related to the data in those rows.
      #lines.import<-readLines(model.stats)[unique(unlist(lines.imp))]
      #params.pi<-names(which(sapply(unique(nameofinterest),grep.fun,main=lines.import)>0))

      #Find the columns that contain these parameters
      index<-list()
      for(l in 1:length(parameters)){
        index[[l]]<-grep(parameters[l],colnames(PP.sample))
      }
      of.interest<-unlist(index)
      #PSA.mat<-as.matrix(PP.sample[,unlist(index)])
      evi<-evppi(of.interest,PP.sample,he)
    }



    ####Calculate EVSI by Quadrature####
    var.prepost<-list()
    start<-Sys.time()

    #foreach(q=1:Q,.packages="R2OpenBugs",.export=ls(.GlobalEnv)) %do% {
    for(q in 1:Q){

      if(exists("ordered",mode = "numeric")){#Creating list of the future data to give to jags
        length.data<-dim(ordered)[2]
        names.data<-data#colnames(ordered)

        Data.Fut<-list()
        for(d in 1:length(index.data)){
          Data.Fut[[d]]<-as.numeric(ordered[q,index.data[[d]]])
        }
        names(Data.Fut)<-names.data
      }
      if(cl.dat=="list"){
        Data.Fut<-data[[q]]
      }

      #Both data sets must be given to jags
      data.full<-append(Data.Fut,data.stats)

      Model.BUGS<- bugs(data.full,inits=NULL,parameters.to.save=moniter,
                        model.file=model.stats, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chains=1,
                        DIC=FALSE,debug=FALSE)
      #Create a dataframe containing all the BUGS samples
      sample<-as.data.frame(Model.BUGS$sims.matrix)

      #Use the BUGS output as inputs for the effects and costs functions
      int<-matrix(NA,nrow=dim(sample)[1],ncol=2*he$n.comparisons)

      #Calculate costs and effects for all simulations
      for(i in 1:dim(sample)[1]){
        #Set function arguements
        formals(effects)<-find.args(effects,sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        #Calculate the incremental effects.
        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,sample,i)
        model.costs<-costs()
        int[i,]<-cbind(model.effects[-he$ref]-model.effects[he$ref],model.costs[-he$ref]-model.costs[he$ref])
      }
      #Calculate the preposterior variance matrix for the costs and effects to calculate EVSI by WTP
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))   }
      print(paste("Update",q,"completed"))
      var.prepost[[q]]<-var(int)
    }
    }

  #Calcualte the EVSI accross different WTP
  mean.var<-apply(simplify2array(var.prepost),1:2,mean)
  EVSI<-array(NA,dim=c(1,length(he$k),1))
  k.<-1
  #EVSI[,1:length(he$k),]<-foreach(k=he$k,.combine="c") %do% {
  for(k in he$k){
    #Fitted INB for each k
    INB<-k*evi$fitted.effects[,1]-evi$fitted.costs[,1]
    #Mean preposterior variance for each k
    Var<-k^2*mean.var[1,1]+mean.var[2,2]-2*k*mean.var[1,2]
    #Rescale fitted INB
    INB.star<-(INB-mean(INB))/sd(INB)*sqrt(max(0,var(he$ib[k.,])-Var))+mean(INB)
    #Calculate EVSI
    EVSI[,k.,]<-mean(pmax(INB.star,0))-max(mean(INB.star),0)
    k.<-k.+1
  }
  #Return EVSI, plus evppi object and bcea object to plot EVSI plus attrib which fits in with later objects..
  to.return<-list(evsi = EVSI,
                  attrib=list(wtp=he$k,N="-",CI="No Uncertainty"),
                  evppi=evi,
                  he=he)
  class(to.return)<-"evsi"
  return(to.return)
}

##comp.evsi.N############################################################
comp.evsi.N<-function(model.stats,data,N,N.range=c(30,1500),effects,costs,he=NULL,parameters=NULL,evi=NULL,Q=50,length.wtp=150,
                      data.stats=NULL,update=c("BUGS","jags"),
                      n.burnin=1000,n.thin=1,n.iter=5000){
  ##'Function generates the data points that will be used to calculate the EVSI for different
  ##'sample sizes.
  ##INPUTS.
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or jags model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. If the future data has already been generated then
  ##'   the data arguement can be given as a list. The elements of the data list should
  ##'   be data lists for jags or BUGS.
  ##'@param N A string defining the name of the variable controlling the sample size in the
  ##'   model.
  ##'@param N.range A two-vector defining the minimum and maximum values considered for the
  ##'   sample size. If the future data is given in the data arguement then N.range should
  ##'   give the sample sizes for which the data was generated.
  ##'@param effects This can either be given as a string which defines the name of the
  ##'   effectivness measure in the BUGS/jags model. Or it can be given as a
  ##'   function that takes the output of the model and calculates the
  ##'   effectiveness measure. The inputs of this function should be parameter
  ##'   names found in the model file.
  ##'@param costs This has a similar format to e but gives the value of the costs in
  ##'   the model file or a function defining the costs.
  ##'@param he A bcea object containing the base case analysis for the model.
  ##'@param parameters A list of the names of the parameters that are being targetted by the study.
  ##'  This can be given OR an evppi object needs to be given.
  ##'@param evi An evppi object that contains the EVPPI analysis for the parameters
  ##'   that are being informed by the sample information.
  ##'@param Q The number of quadrature points used the estimate the EVSI.
  ##'@param length.wtp The number of willingness to pay values that should be considered to estimate the EVSI
  ##'@param data.stats A data file for the BUGS/jags model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the jags/BUGS model
  ##'@param n.burnin The burnin for the jags/BUGS model
  ##'@param n.iter The number of interations for the jags/BUGS model
  ##'
  ##'OUTPUTS.
  ##'@return An evsi.comp object to be used in the evsi.calc function.
  ##'...


  #Is data in the correct format??
  cl.dat<-class(data)
  if(cl.dat=="list"){
    if(length(data)!=length(N.range)){stop("The number of future data samples is not equal to the number of sample sizes.")}

    Q<-length(N.range)
    if(Q<30){warning("You are calculating the EVSI with fewer than 30 future samples. The EVSI is unlikely to be accurately estimated.")}
  }

  #Set the functions for costs and effects
  if(class(costs)!="function"){
    costs<-function(costs){
      return(costs)
    }
  }

  if(class(effects)!="function"){
    effects<-function(effects){
      return(effects)
    }
  }

  #Extract the arguements of the effects and costs functions
  nameofinterest<-unique(c(names(formals(effects)),names(formals(costs))))

  #Check which of these parameters are also in the model file
  #grep reads the model file finds the occurance of a string.
  grep.fun<-function(name,main){
    length(grep(name, main, value = FALSE))
  }

  #Set the function arguements for the model functions
  find.args<-function(func,sample,i){
    #Find the names of all the function arguements
    args.names<-names(formals(func))
    #Find the number of times a variable is containted in the jags model
    #This allows us to determine matrix contructions for certain variables.
    #The variable a vector - i.e. it only appears once in the jags output
    jags.params<-which(sapply(moniter,grep.fun,main=names(sample))==1)
    inputs<-sample[,moniter[jags.params]]

    #Matrix Parameters
    multiple.params<-moniter[-jags.params]

    #Function to extract the matrix parameters from the jags object
    multiple.params.extract<-function(multi.params){
      if(length(multi.params)==0){return(NULL)}
      else{
        list.params<-list()
        length(list.params)<-length(multi.params)
        for(j in 1:length(multiple.params)){
          list.params[j]<-list(as.numeric(sample[i,grep(multi.params[j],names(sample))]))
        }
        return(list.params)}
    }

    #Sets the effects function arguements for the single parameters
    formals(func)[which(args.names %in% names(inputs))]<-
      inputs[i,which(names(inputs) %in% args.names)]
    #Sets the effects function arguements for the matrix parameters
    formals(func)[args.names %in% multiple.params]<-
      multiple.params.extract(args.names[which(args.names %in% multiple.params)])
    return(formals(func))
  }

  quantiles<-sample(seq(0,1,length.out=Q),replace=F)
  if(length(N.range)==2){N.samp<-trunc((seq(sqrt(N.range[1]),sqrt(N.range[2]),length.out = Q))^2)}
  if(length(N.range)>2){N.samp<-N.range}

  #Generate future samples by finding prior-predictive distribution
  if(update=="jags"){#Model can be written in jags.
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    ####Calculate EVSI by Quadrature####
    var.prepost<-list()
    start<-Sys.time()
    for(q in 1:Q){
      #Set the data including the sample size.
      Samp.Size<-list(N.samp[q])
      names(Samp.Size)<-N

      if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
        #Runs the jags model based on the rjags package
        Model.JAGS<- jags.model(model.stats,data=append(Samp.Size,data.stats),quiet=TRUE)
        update(Model.JAGS,n.burnin,progress.bar="none")
        Prior.Pred <- coda.samples(Model.JAGS, prior.pred.data, n.iter=n.iter,n.thin=n.thin,progress.bar="none")
        PP.sample<-as.data.frame(Prior.Pred[[1]])
      }

      if(cl.dat=="character"){
        #Determine which columns contain the data
        index.data<-list()
        for(l in 1:length(data)){
          index.data[[l]]<-grep(data[l],colnames(PP.sample))
        }
        length.data<-length(unlist(index.data))
        names.data<-data
        Data.Fut<-array()
        for(d in 1:length.data){
          #Creating list of the future data to give to jags
          Data.Fut[d]<-as.numeric(quantile(PP.sample[,unlist(index.data)[d]],probs=quantiles[q]))
        }
        Data.Fut.list<-list()
        for(d in 1:length(index.data)){
          Data.Fut.list[[d]]<-Data.Fut[index.data[[d]]]
        }
        Data.Fut<-Data.Fut.list
        names(Data.Fut)<-names.data
      }
      if(cl.dat=="list"){
        Data.Fut<-data[[q]]
      }

      if(class(he)!="bcea"){
        #Calculate costs and effects for all simulations
        i<-1
        #Set function arguements
        formals(effects)<-find.args(effects,PP.sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        e<-c<-matrix(NA,nrow=dim(PP.sample)[1],ncol=length(model.effects))
        for(i in 1:dim(PP.sample)[1]){
          #Set function arguements
          formals(effects)<-find.args(effects,PP.sample,i)
          #Runs model with all arguements set as above
          e[i,]<-effects()
          #Repeat analysis for costs
          #Set function arguements
          formals(costs)<-find.args(costs,PP.sample,i)
          c[i,]<-costs()
        }
        he<-bcea(e,c)
      }
      #Check number of interventions
      if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
                                  multi-decision problems")}

      if(class(evi)!="evppi"){
        #Determine which model rows contain statements defining the data
        #lines.imp<-list()
        #for(l in 1:length(data)){
        #  lines.imp[[l]]<-grep(data[l],readLines(model.stats))
        #}

        #Determine which parameters are related to the data in those rows.
        #lines.import<-readLines(model.stats)[unique(unlist(lines.imp))]
        #params.pi<-names(which(sapply(unique(nameofinterest),grep.fun,main=lines.import)>0))

        #Find the columns that contain these parameters
        index<-list()
        for(l in 1:length(parameters)){
          index[[l]]<-grep(parameters[l],colnames(PP.sample))
        }
        of.interest<-unlist(index)
        #PSA.mat<-as.matrix(PP.sample[,unlist(index)])
        evi<-evppi(of.interest,PP.sample,he)
      }



      #Both data sets must be given to jags
      data.full<-append(append(Data.Fut,data.stats),Samp.Size)

      Model.JAGS<- jags.model(model.stats,data =  data.full,quiet = TRUE)
      update(Model.JAGS,n.burnin,progress.bar="none")
      samples <- coda.samples(Model.JAGS, moniter, n.iter=n.iter,n.thin=n.thin,n.chain=1,progress.bar="none")
      #Create a dataframe containing all the JAGS samples
      sample<-as.data.frame(samples[[1]])

      #Use the JAGS output as inputs for the effects and costs functions
      e.int<-c.int<-array()
      length(e.int)<-length(c.int)<-dim(sample)[1]

      #Calculate costs and effects for all simulations
      for(i in 1:dim(sample)[1]){
        #Set function arguements
        formals(effects)<-find.args(effects,sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        #Calculate the incremental effects.
        e.int[i]<-model.effects[-he$ref]-model.effects[he$ref]

        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,sample,i)
        model.costs<-costs()
        c.int[i]<-model.costs[-he$ref]-model.costs[he$ref]
      }
      #Calculate the preposterior variance matrix for the costs and effects to calculate EVSI by WTP
      var.prepost[[q]]<-var(cbind(e.int,c.int))
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))
      }
      print(paste("Update",q,"completed"))
      }
  }
  #Analysis in BUGS
  if(update=="BUGS"){#Model can be written in BUGS.
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    var.prepost<-list()
    start<-Sys.time()
    for(q in 1:Q){
      #Set the data including the sample size.
      Samp.Size<-list(N.samp[q])
      names(Samp.Size)<-N

      if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
        #Runs the BUGS model based on the R2OpenBUGS package
        Model.BUGS<- bugs(data.stats,inits=NULL,parameters.to.save=append(prior.pred.data,Samp.Size),
                          model.file=model.stats, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chain=1,
                          DIC=FALSE,debug=FALSE)
        PP.sample<-as.data.frame(Model.BUGS$sims.matrix)
      }


      if(cl.dat=="character"){
        #Determine which columns contain the data
        index.data<-list()
        for(l in 1:length(data)){
          index.data[[l]]<-grep(data[l],colnames(PP.sample))
        }
        length.data<-length(unlist(index.data))
        names.data<-data
        Data.Fut<-array()
        for(d in 1:length.data){
          #Creating list of the future data to give to jags
          Data.Fut[d]<-as.numeric(quantile(PP.sample[,unlist(index.data)[d]],probs=quantiles[q]))
        }
        Data.Fut.list<-list()
        for(d in 1:length(index.data)){
          Data.Fut.list[[d]]<-Data.Fut[index.data[[d]]]
        }
        Data.Fut<-Data.Fut.list
        names(Data.Fut)<-names.data
      }
      if(cl.dat=="list"){
        Data.Fut<-data[[q]]
      }

      if(class(he)!="bcea"){
        #Calculate costs and effects for all simulations
        i<-1
        #Set function arguements
        formals(effects)<-find.args(effects,PP.sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        e<-c<-matrix(NA,nrow=dim(PP.sample)[1],ncol=length(model.effects))
        for(i in 1:dim(PP.sample)[1]){
          #Set function arguements
          formals(effects)<-find.args(effects,PP.sample,i)
          #Runs model with all arguements set as above
          e[i,]<-effects()
          #Repeat analysis for costs
          #Set function arguements
          formals(costs)<-find.args(costs,PP.sample,i)
          c[i,]<-costs()
        }
        he<-bcea(e,c)
      }
      #Check number of interventions
      if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
                                  multi-decision problems")}

      if(class(evi)!="evppi"){
        #Determine which model rows contain statements defining the data
        #lines.imp<-list()
        #for(l in 1:length(data)){
        #  lines.imp[[l]]<-grep(data[l],readLines(model.stats))
        #}

        #Determine which parameters are related to the data in those rows.
        #lines.import<-readLines(model.stats)[unique(unlist(lines.imp))]
        #params.pi<-names(which(sapply(unique(nameofinterest),grep.fun,main=lines.import)>0))

        #Find the columns that contain these parameters
        index<-list()
        for(l in 1:length(parameters)){
          index[[l]]<-grep(parameters[l],colnames(PP.sample))
        }
        of.interest<-unlist(index)
        #PSA.mat<-as.matrix(PP.sample[,unlist(index)])
        evi<-evppi(of.interest,PP.sample,he)
      }

      #Both data sets must be given to jags
      data.full<-append(append(Data.Fut,data.stats),Samp.Size)
      ####Calculate EVSI by Quadrature####

      Model.BUGS<- bugs(data.full,inits=NULL,parameters.to.save=moniter,
                        model.file=model.stats, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chains=1,
                        DIC=FALSE,debug=FALSE)
      #Create a dataframe containing all the BUGS samples
      sample<-as.data.frame(Model.BUGS$sims.matrix)

      #Use the BUGS output as inputs for the effects and costs functions
      e.int<-c.int<-array()
      length(e.int)<-length(c.int)<-dim(sample)[1]

      #Calculate costs and effects for all simulations
      for(i in 1:dim(sample)[1]){
        #Set function arguements
        formals(effects)<-find.args(effects,sample,i)
        #Runs model with all arguements set as above
        model.effects<-effects()
        #Calculate the incremental effects.
        e.int[i]<-model.effects[-he$ref]-model.effects[he$ref]

        #Repeat analysis for costs
        #Set function arguements
        formals(costs)<-find.args(costs,sample,i)
        model.costs<-costs()
        c.int[i]<-model.costs[-he$ref]-model.costs[he$ref]
      }
      #Calculate the preposterior variance matrix for the costs and effects to calculate EVSI by WTP
      var.prepost[[q]]<-var(cbind(e.int,c.int))
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))
      }
      print(paste("Update",q,"completed"))
      } }

  #Calcualte the EVSI accross different WTP
  print("Using curve fitting to find EVSI for alternative sample sizes")
  wtp.seq<-seq(min(he$k),max(he$k),length.out=length.wtp)
  model.ab<-function(){
    beta~dnorm(Nmax/2,shape.Nmax)%_%T(0,)
    for(i in 1:N){
      y[i]~dnorm(mu[i],tau)
      mu[i]<-var.PI*(x[i]/(x[i]+beta))
    }
    sigma~dt(sigma.mu,sigma.tau,3)%_%T(0,)#dunif(0.01,50)
    tau<-1/sigma^2
  }
  file.curve.fitting <- file.path("~",fileext="model_curve_fitting_EVSI.txt")
  write.model(model.ab,file.curve.fitting)

  beta.mat<-array(NA,dim=c(3000,length(wtp.seq)))
  k.<-1
  start<-Sys.time()
  for(k in wtp.seq){
    Var.X.prob<-k^2*simplify2array(var.prepost)[1,1,]+simplify2array(var.prepost)[2,2,]-
      2*k*simplify2array(var.prepost)[1,2,]

    fitted.wtp<-k*evi$fitted.effects[,1]-evi$fitted.costs[,1]
    full.wtp<-k*he$delta.e-he$delta.c
    y<-t(var(full.wtp)-Var.X.prob)

    data.a.b<-list(sigma.mu=sd(y)/2,
                   sigma.tau=1/sd(y),
                   N=length(N.samp),
                   shape.Nmax=0.0005/max(N.samp),
                   var.PI=var(fitted.wtp),
                   Nmax=max(N.samp),
                   y=as.vector(y),
                   x=as.vector(N.samp)
    )

    n.burnin <- 1000  # Number of burn in iterations
    n.thin<-1
    n.iter <- 3000 # Number of iterations per chain

    # Perform the MCMC simulation with JAGS.
    #Runs the jags model based on the rjags package
    Model.JAGS<- jags.model(file.curve.fitting,data=data.a.b,quiet=TRUE)
    update(Model.JAGS,n.burnin,progress.bar="none")
    beta.ab <- coda.samples(Model.JAGS, c("beta"), n.iter=n.iter,n.thin=n.thin,progress.bar="none")

    beta.mat[,k.]<-as.data.frame(beta.ab[[1]])[,1]

    if(k.==1){end<-Sys.time()
    comp.time<-end-start
    print(paste(c("Curve fitting updating requires ",round(comp.time,2)," seconds. There are ",length.wtp," willingness-to-pay values to consider.
                  The remaining computation time is around ",round(comp.time*(length.wtp-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                sep="",collapse = ""))
    }
    print(paste("Update",k.,"completed"))
    k.<-k.+1
  }
  #Return EVSI, plus evppi object and bcea object to plot EVSI.
  to.return<-list(beta = beta.mat,
                  N=N.samp,
                  wtp=wtp.seq,
                  evi=evi,
                  he=he)
  class(to.return)<-"evsi.N"
  return(to.return)

  }

##evsi.calc############################################################
evsi.calc<-function(comp.evsi.N,wtp=NULL,N=NULL,CI=NULL){
  ##'Calculates the EVSI accross willingness-to-pay, sample size and gives uncertainty for these
  ##'measurements as required. This function can be fed into/used in the other plot functions
  ##'This should allow people to upload their own EVSI matrices and still use the plotting capabilities.
  ##INPUTS
  ##'@param comp.evsi.N Output of the comp.evsi.N function
  ##'@param wtp The willingness to pay value(s) for which the EVSI should be calculated.
  ##'   It will be chosen as all the wtp values if NULL.
  ##'@param N The sample sizes for which the EVSI should be calculated. If NULL it will select
  ##'   all the values that the EVSI has been calculated for. NOTE: N can be chosen as any
  ##'   value but wtp needs to be on the grid chosen in comp.evsi.N.
  ##'@param CI The confidence levels for the credible intervals of the EVSI. If NULL chosen as
  ##'   c(0.025,0.25,0.5,0.75,0.975)
  ##'
  ##'OUTPUTS
  ##'@return An evsi object.
  ##'1. evsi An array containing the EVSI by wtp, N and across different uncertaincies
  ##'2. attrib A list of wtp, N and prob describing the attributes of the evsi matrix.
  ##'3. evppi An evppi object containing all the information about the calculation of
  ##'   EVPPI.
  ##'4. he A bcea object containing all the information about the underlying health
  ##'   economic model
  ##'   @example
  ##'   ...

  if(class(comp.evsi.N)!="evsi.N"){stop("comp.evsi.N must be calculated using the comp.evsi.N function. Please use this function. evsi objects can be created using the evsi.upload function.")}
  #Format wtp
  cl<-class(wtp)
  #Select all wtp values if NULL
  if(cl!="numeric"){
    wtp<-comp.evsi.N$wtp
  }
  wtp.length<-length(wtp)
  #Select wtp from the grid of wtp if chosen by user.
  if(cl=="numeric"){
    for(i in 1:wtp.length){
      wtp[i]<-comp.evsi.N$wtp[which.min(abs(wtp[i]-comp.evsi.N$wtp))]
    }
  }

  #Format N
  cl<-class(N)
  #Select all N values if NULL
  if(is.null(N)){
    N<-seq(min(comp.evsi.N$N),max(comp.evsi.N$N),by=10)
  }
  N.length<-length(N)

  #Format CI
  cl<-class(CI)
  #Set the default CI values if NULL
  if(cl!="numeric"){
    CI<-c(0.025,0.25,0.5,0.75,0.975)
  }
  CI.length<-length(CI)
  #Set N as requested by user
  if(cl=="numeric"){
    CI<-CI
  }
  #Use CI to create a beta matrix
  #Find the appropriate quantiles for the beta parameter of interest
  beta.quantiles<-apply(comp.evsi.N$beta,2,quantile,prob=CI)

  #Extracting the beta parameter for the wtp of interest
  beta.focal<-as.matrix(beta.quantiles[,which(comp.evsi.N$wtp%in%wtp)])

  EVSI.mat<-array(NA,dim=c(N.length,wtp.length,CI.length))

  non.zero<-which(colSums(comp.evsi.N$evi$fitted.effects)!=0)
  calc.EVSI<-function(beta.focal,wtp,N){
    #Find the fitted values for the wtp
    fitted.PI<-wtp*comp.evsi.N$evi$fitted.effects[,non.zero]-comp.evsi.N$evi$fitted.costs[,non.zero]
    #Variance for a specific N and beta
    var.pre<-var(fitted.PI)*(N/(N+beta.focal))
    #Rescaled fitted values
    fitted.star<-(fitted.PI-mean(fitted.PI))/sd(fitted.PI)*sqrt(var.pre)+mean(fitted.PI)
    #Calculate EVSI
    EVSI<-mean(pmax(fitted.star,0))-max(mean(fitted.star),0)
    return(EVSI)
  }

  for(i in 1:CI.length){
    for(j in 1:wtp.length){
      EVSI.mat[,j,i]<-sapply(N,calc.EVSI,beta.focal=beta.focal[i,j],wtp=wtp[j])
    }
  }
  to.return<-list(evsi=EVSI.mat,
                  attrib=list(wtp=wtp,N=N,CI=CI),
                  evppi=comp.evsi.N$evi,
                  he=comp.evsi.N$he)
  class(to.return)<-"evsi"
  return(to.return)
}

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
  he<-bcea(e,c)

  #Calculate EVPPI object
  evi<-evppi(parameter,input,he,method="gam")

  to.return<-list(evsi=EVSI.arr,attrib=list(wtp=wtp,N=N,CI=0.5),evppi=evi,he=he)
  class(to.return)<-"evsi"
  return(to.return)
}


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
                                cex = 0.7, bty = "n", lty = c(1, 1,1),
                                lwd = c(2, 1))}
  if(select.length>1){legend(alt.legend,legend=c(min(N),rep(NA,max(0,select.length-2)),max(N)),
                             fill=colours,border=colours,cex=0.75,
                             y.intersp=max(0.1,1.2/select.length),bty="n")}
}

##evsi.pop############################################################
evsi.pop<-function(evsi,trial.cost=NULL,setup=NULL,pp=NULL,
                   Pop=c(0,10000),Time=c(1,20),Dis=0.035,
                   wtp=NULL,N=NULL,pos=c("topright")){
  ##'Produces a graphic that considers the cost-effectiveness of a single future trial for different
  ##'Time Horizons for the treatment and incidence population
  ##'INPUTS:
  ##'@param evsi An evsi object.
  ##'@param trial.cost Either a fixed trial cost (numeric) or a distribution for a trial cost (character string)
  ##'   - to give probability of cost-effectiveness for the trial. Can be NULL if setup and pp are both given
  ##'@param setup Gives the setup costs of the trial. Can be NULL if trial.costs are given. If NULL and pp given
  ##'   setup is taken to be 0.
  ##'@param pp Gives the per person costs of the trial. Can be NULL if trial.costs are given.
  ##'@param Pop The minimum and maximum possible values for the number of people benefitting from
  ##'   the treatment
  ##'@param Time The minimum and maximum possible values for the time horizen of the treatment.
  ##'@param Dis The discount rate for future studies, based on NICE guidelines.
  ##'@param wtp The willingness to pay value that this analysis is being undetaken for.
  ##'   If NULL then the function will automatically select the central wtp to default in BCEA will be 25000
  ##'@param N The sample size which we want to consider for the trial. If NULL then taken as the median trial size
  ##'   considered.
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
  if((class(N)!="numeric")||(class(evsi$attrib$N)!="numeric")){
    N.select<-ceiling(length(evsi$attrib$N)/2)
    N<-evsi$attrib$N[N.select]
  }
  if(class(N)=="numeric"){
    N.select<-which.min((evsi$attrib$N-N)^2)
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
    evsi.1<-evsi$evsi[N.select,wtp.select,length(evsi$attrib$CI)]
    evsi.2<-evsi$evsi[N.select,wtp.select,length(evsi$attrib$CI)-1]
    evsi.params<-c((q.2*evsi.1-evsi.2*q.1)/(q.2-q.1),(evsi.2-evsi.1)/(q.2-q.1))
    quants<-array()
    for(i in 1:(length(evsi$attrib$CI)-1)){
      quants[i]<-pnorm(evsi$evsi[N.select,wtp.select,i],evsi.params[1],evsi.params[2])-evsi$attrib$CI[length(evsi$attrib$CI)-(i-1)]}
    discrep<-mean(quants,na.rm=T)
    if(discrep>0.2){warning("EVSI distribution not well approximated by a normal, the probabilistic EVSI analysis may be incorrect")}
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
    if(length(trial.cost)==1){trial.cost<-c(trial.cost,0)}

    Time.min<-min(Time);Time.max<-max(Time);Pop.min<-min(Pop);Pop.max<-max(Pop)
    dens.points<-100
    Time.seq<-seq(Time.min,Time.max,length.out=dens.points)
    Pop.seq<-seq(Pop.min,Pop.max,length.out=dens.points)

    Prob.mat<-matrix(NA,nrow=length(Time.seq),ncol=length(Pop.seq))
    for(i in 1:length(Time.seq)){
      for(j in 1:length(Pop.seq)){
        ENBS.mean<-ENBS(evsi.params[1],Pop.seq[j],Time.seq[i],Dis,trial.cost[1])
        ENBS.sd<-ENBS.sd.calc(evsi.params[2],Pop.seq[j],Time.seq[i],Dis,trial.cost[2])
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
        main="Probability of Cost-Effective Trial",xlab="Time Horizen",ylab="Incidence Population",
        xlim = c(Time.min,Time.max),ylim=c(Pop.min,Pop.max),
        breaks=seq(0,1,length.out=101))

  legend(alt.legend,c("Prob=0",rep(NA,98/2),"Prob=.5",rep(NA,96/2),"Prob=1"),fill=colours,border=colours,cex=0.75,y.intersp=0.15)

}

##plot.evsi.N###########################################################
plot.evsi.N<-function(evsi,wtp=NULL,pos=c("bottomright"),CI=NULL){
  ##'Calculating the EVSI for a specific WTP giving the uncertainty bands across the different
  ##'samples sizes
  ##INPUTS
  ##'@param comp.evsi.N Output of the comp.evsi.N function
  ##'@param wtp The willingness to pay value that the graphic should be produced for - it will
  ##'   be chosen if wtp=NULL.
  ##'@param prob The confidence levels for the uncertainty bands.
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
       col="white",xlab=expression("Sample Size"),ylab="Per Person EVSI",oma=c(0,0,-1,0))

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

  fitted.PI<-wtp*evsi$evppi$fitted.e-evsi$evppi$fitted.c
  abline(h=mean(apply(fitted.PI,1,max))-max(apply(fitted.PI,2,mean)),col="springgreen",lwd=lwd[CI.length+1],lty=lty[CI.length+1])

  legend(alt.legend,c(as.character(CI),"EVPPI"),
         col=c(rep("black",CI.length),"springgreen"),lwd=lwd,lty=lty)

}

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
  ##OUTPUTS
  ##'@return SS The optimal sample size.
  ##'@return SS.CI Credible intervals for the optimal sample size. TO BE IMPLEMENTED

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
  if((max.less<1)||(max.greater>length(ENBS))){
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

  return(list(SS.max=N.max,ENBS=ENBS.max))

}

##evsi.plot.dynam######################################################
evsi.plot.dynam<-function(evsi,width=600,height=600){
  ##'Function launches shiny app to show dynamic graphics
  ##'INPUTS
  ##'@param evsi An evsi object for which dynamic graphics are needed.
  ##'
  ##OUTPUTS
  ##'@return A dynamic shiny app that plots the EVSI for different wtp,
  ##'sample size, incidence population, time horizon and costs as required.
  digits<-nchar(as.character(round(diff(evsi$attrib$wtp))))[1]
  if(is.na(digits)){rounding<-0}
  else{rounding<-1-digits}
  ui<-shiny::fluidPage(theme=shinytheme('united'),
                       titlePanel(h1("Visualisations for the EVSI")),
                       mainPanel(
                         #Create Tabs
                         tabsetPanel(
                           #Description Panel
                           tabPanel("Introduction to plotting",
                                    fluidRow(p("The Expected Value of Sample Information can be used to determine the economic benefit of a future
                                               trial. However, in general, this economic benefit is conditional on several deterministic inputs that can
                                               generally have a large impact on the value of the EVSI and any decision you make with the EVSI in mind.
                                               Therefore, we have developed this tool to allow for the simple visualistion of the EVSI along with easy
                                               manipulation of these deterministic inputs. This tool also allows for a simple user interface that can be
                                               shown to collaborators and other stakeholders.",
                                               h2("Willingness to Pay"),
                                               p("The willingness to pay (WTP) is the amount of money the decision maker has avaliable to pay for an
                                                 additional 1 unit of benefit. This is typically given as a range and therefore the EVSI can be visualised
                                                 across different WTP values in the tab WTP. For the other analyses the WTP must be fixed. Although the plots
                                                 can easily be redrawn for different values of the WTP allowing the analysis to be easily completed for a
                                                 large number of different thresholds."),
                                               h2("Sample Size"),
                                               p("The sample size of the future trial is rarely determined before the EVSI analysis and one of the analysis
                                                 avaliable in this tool involves determining the optimal sample size for the future trial. The EVSI can be
                                                 visualised by sample size on the tab N. The value of a sample is bounded above by the value of resolving
                                                 all uncertainty in the parameters of interest. The key information to be gleaned from the by N plot is the
                                                 speed at which the EVSI reachs this upper bound."),
                                               h2("Trial Cost-effectivenes"),
                                               p("A key use of the EVSI is to determine whether a future trial will be cost-effective. This means that the
                                                 value of the future trial exceeds the cost of the trial. While in general this is a simple extension of the
                                                 ideas underpinning the EVSI, it typically depends on some additional inputs which are rarely known with
                                                 certainty."))),
                                    fluidRow(column(11,offset=1,h3("Trial costs"),
                                                    p("The cost of the trial must be specified to allow comparison with the EVSI. However, these are rarely known
                                                      with certainty and so this application allows a range of values to be specified for the cost and then
                                                      considers the cost-effectiveness of the trial taking into account this uncertainty."),
                                                    h3("Incidence Population"),
                                                    p("The EVSI is calculated as the EVSI per person who will benefit from the treatments under consideration.
                                                      Therefore, to compare with the trial costs, we must multiply by the number of people who will benefit
                                                      from the treatments. We call this the \"Incidence Population\". While there may be some literature that
                                                      can inform this parameter it cannot be known and therefore the cost-effectiveness of the trial is determined
                                                      for different levels of this population. This is an important consideration as the level of this population
                                                      can make a large difference to the cost-effectiveness of the trial."),
                                                    h3("Time Horizon"),
                                                    p("Finally, the time horizon of a treatment is the length of time the treatments will be available in the
                                                      market before an alternative superseeds the treatments by being more cost-effective. Typically, this is
                                                      assumed to be around 10 years but we allow variation in this parameter to consider the cost-effectiveness
                                                      for alternative time horizons. This is powerful as it is rarely known the length of time a technology will
                                                      be available in the market.")))
                                                    ),

                           #Per person EVSI by Willingness to Pay - n slider
                           tabPanel("EVSI by Willingness To Pay",
                                    sidebarPanel(p("The EVSI changes depending on the willingness-to-pay of the underlying decision maker and therefore
                                                   it is useful to visualise the EVSI for different WTP thresholds. This plot is also used visualise the
                                                   EVPI - which gives the maximum value for ANY future trial - and the EVPPI - which gives a maximum
                                                   value for a study targeting uncertainty in the parameters of interest of our study.
                                                   The EVSI also changes for different sample sizes as studies increase in value as the sample size increases.
                                                   The sample size can therefore be changed:"),
                                                 selectInput(inputId="n",label="Choose a sample size",
                                                             choices=evsi$attrib$N,
                                                             selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)]),
                                                 p("The EVSI should remain below the EVPPI for all sample sizes and willingness-to-pay thresholds.
                                                   All three measures should reach a sharp peak - this represents the \"break-even\" point between the
                                                   two treatment options where the decision uncertainty is at its maximum as the two treatment options are
                                                   equally likely to be cost-effective."),width=4),
                                    mainPanel(plotOutput(outputId="ppEVSIbywtp",width=width,height=height),
                                              p("Citations: McCabe, C., Claxton K. and Culyer, A., ",
                                                a("The NICE cost-effectiveness threshold",href="https://link.springer.com/article/10.2165/00019053-200826090-00004"),
                                                em(", PharmacoEconomics"),", 2008"),
                                              p("Baio, G., ",
                                                a("Bayesian Methods in Health Economics",href="https://sites.google.com/a/statistica.it/gianluca/bookhe"),
                                                em(", Springer"),", 2012"),
                                              p("Baio, G., Berardi, A. and Heath, A., ",
                                                a("Bayesian Cost-Effectiveness Analysis with the R package BCEA",href="http://www.springer.com/gb/book/9783319557168"),
                                                em(", Springer"),", 2017"),width=8)
                                                 ),

                           #Per person EVSI by N - wtp slider

                           tabPanel("EVSI by Sample Size",
                                    sidebarPanel(p("The EVSI increases as the sample size of the underlying trial increases. This graphic shows the
                                                   EVSI across different sample sizes. This relationship with N changes depending on the value of the
                                                   willingness-to-pay so the plot can be considered for changing values of the WTP.
                                                   The EVSI is calculated using Bayesian regression. This means that posterior credible intervals for the EVSI
                                                   can be calculated and these are plotted on the graphic. This demonstrates the uncertainty in the EVSI
                                                   estimate. If the uncertainty is too large for decision making then Q should be increased in the
                                                   EVSI calculation."),
                                                 selectInput(inputId="wtp",label="Choose a Willingness-to-Pay Threshold",
                                                             choices=round(evsi$attrib$wtp,rounding),
                                                             selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar)^2)],rounding)),
                                                 #sliderInput(inputId="wtp",label="Choose a Willingness-to-Pay Threshold",
                                                 #min=min(evsi$attrib$wtp),max=max(evsi$attrib$wtp),value=min(evsi$he$kstar),
                                                 #step=round(evsi$attrib$wtp[2]-evsi$attrib$wtp[1])),
                                                 p(paste("In general, the EVSI will be more accurately estimated for higher values of the EVSI. This
                                                         would will be for willingness to pay values close to the \"break-even\" point of ",
                                                         round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar)^2)],rounding),".")),width=4),
                                    mainPanel(plotOutput(outputId="ppEVSIbyn",width=width,height=height),
                                              p("Citations: Heath, A., Manolopoulou I. and Baio, G., ",
                                                a("Efficient Monte Carlo Estimation of the Expected Value of Sample Information using Moment Matching",href="https://arxiv.org/abs/1611.01373"),
                                                em(", arXiv Preprint"),", 2017"),
                                              p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                em(", Working Paper"),", 2017"),
                                              width=8)
                                    ),

                           #Trial Costs Input
                           tabPanel("Cost-effectiveness of a Trial",column(12,
                                                                           tabsetPanel(tabPanel("Trial Costs",fluidRow(
                                                                             column(5,
                                                                                    fluidRow(column(12,p("To determine the cost-effectiveness of a trial, the EVSI must be compared with the costs of
                                                                                                         undertaking the trial. In general, these costs are split into two categories; the setup and the per person costs.
                                                                                                         The are setup costs are the overhead costs of the trial and will be incurred irrespective of the size
                                                                                                         of the trial. Typical setup costs may be training for staff or the purchase of specialised equipment.
                                                                                                         Additional costs will then be incurred by each participant in the trial such as the cost of
                                                                                                         administering the treatment or following up the patient."),
                                                                                                    p("It is unlikely that these costs will be known exactly. Therefore, you should give a range of possible
                                                                                                      values for the costs. If the costs are known with certainty then simply input the same values for
                                                                                                      the maximum and minimum possible values of the costs. ")
                                                                                                    )),
                                                                                    fluidRow(
                                                                                      column(6,numericInput(inputId="Setupmin",label="Minimum Setup Costs for the Trial",
                                                                                                            value=1000,step=10,min=0),
                                                                                             numericInput(inputId="PerPersmin",label="Minimum Cost Per Person",min=0,
                                                                                                          value=1500,step=10)),
                                                                                      column(6,numericInput(inputId="Setupmax",label="Maximum Setup Costs for the Trial",
                                                                                                            value=100,step=10,min=0),
                                                                                             numericInput(inputId="PerPersmax",label="Maximum Cost Per Person",min=0,
                                                                                                          value=1000,step=10)
                                                                                      ))),
                                                                             column(5,fluidRow(p("The cost-effectiveness of a trial depends on two additional inputs. These are the incidence population, i.e.
                                                                                                 the yearly incidence of the disease under consideration. This gives the number of patients that will benefit from
                                                                                                 the treatment in each year that the treatment is used. This may not be known with certainty and the plot considers
                                                                                                 the cost-effectiveness of the trial for different possible numbers of patients. However, it is necessary to give
                                                                                                 possible values of the incidence population."),
                                                                                               p("The time horizon gives the number of years that the most cost-effective treatment will be available.
                                                                                                 This can be thought of as the number of years before a more effective treatment will be developed.
                                                                                                 This will depend on the disease areas as fast moving diseases such as cancer will have a shorter time
                                                                                                 horizon. The maximum and minimum possible values for the time horizon should be specified here.")),
                                                                                    fluidRow(column(6,numericInput(inputId="Popmin",label="Minimum Incidence Population",value=0,step=100,min=0),
                                                                                                    numericInput(inputId="Timemin",label="Minimum Time Horizon",value=0,step=1,min=0)),
                                                                                             column(6,
                                                                                                    numericInput(inputId="Popmax",label="Maximum Incidence Population",value=1e+05,step=100,min=0),
                                                                                                    numericInput(inputId="Timemax",label="Maximum Time Horizon",min=0,value=25,step=1)
                                                                                             ))),
                                                                             column(2,p("The final input to determine the cost-effectiveness of the trial is the discount rate for future treatments.
                                                                                        In general, health benefits now are more valuable than health benefits in the future. NICE recommend 3.5% as the
                                                                                        discount rate for the treatments but this can be changed here."),
                                                                                    numericInput(inputId="Dis",label="Discount Rate",
                                                                                                 value=0.035,step=0.001))),
                                                                             fluidRow(column(12,p("Citations:")))),
                                                                             tabPanel("Probability of CE Trial",
                                                                                      fluidRow(sidebarPanel(#Population EVSI - prob of CE plot
                                                                                        selectInput(inputId="n.CE",label="Choose a sample size",
                                                                                                    choices=evsi$attrib$N,
                                                                                                    selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)]),
                                                                                        selectInput(inputId="wtp.CE",label="Choose a Willingness-to-Pay Threshold",
                                                                                                    choices=round(evsi$attrib$wtp,rounding),selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar)^2)],rounding)),
                                                                                        uiOutput("PopDynam"),
                                                                                        uiOutput("TimeDynam"),width=4),
                                                                                        mainPanel(plotOutput(outputId="ProbCE",width=width,height=height),
                                                                                                  width=8)),
                                                                                      fluidRow(column(12,p("The probability of having a cost effective trial is shown for different time horizons and incidence
                                                                                                           populations. The plot is",em("white"),"when the trial is cost-effective and",em("black"),"when the trial
                                                                                                           is not cost-effective. A trial is more likely to be cost effective as the time horizon and incidence
                                                                                                           population increases and therefore the plot will typically be white in the top right hand corner."),
                                                                                                      p("A trial is",em("cost-effective"),"when the expected value of the information gained in the trial
                                                                                                        is greater than the trial costs. As both the trial costs and the EVSI are not known with certainty
                                                                                                        there are some values for the time horizon and incidence population where the costs of the trial could
                                                                                                        exceed the value but it is uncertain whether this will be the case. In these settings, the plot is
                                                                                                        blue. The darker the blue, the more likely it is that the costs of the trial exceed the EVSI."),
                                                                                                      p("Citations: Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                        a("EVSI Visualisations for cost-effective trial analysis",href="https://arxiv.org/abs/1611.01373"),
                                                                                                        em(", NOT WRITTEN"),", 2018"),
                                                                                                      p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                        a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                                                                        em(", Working Paper"),", 2017")))
                                                                                      ),
                                                                             #Optimal Sample Size
                                                                             tabPanel("Optimal Sample Size",
                                                                                      sidebarPanel(selectInput(inputId="wtp.OS",label="Choose a Willingness-to-Pay Threshold",
                                                                                                               choices=round(evsi$attrib$wtp,rounding),
                                                                                                               selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar)^2)],rounding)),
                                                                                                   uiOutput("Pop.OSDynam"),
                                                                                                   uiOutput("Time.OSDynam"),width=4),

                                                                                      mainPanel(fluidRow(p("It is possible to find the sample size for your trial that will give the maximum value for money.
                                                                                                           While it is possible to find this value for all trials, there is no guarantee that this optimal sample
                                                                                                           size will yield a positive benefit from the sampling. In these cases, the EVSI indicates that the trial
                                                                                                           is not cost-effective and the ",em("Expected Net Benefit of Sampling"),"(ENBS) is less than 0."),
                                                                                                         p("It is important to note that the optimal sample size can only be found between",
                                                                                                           min(evsi$attrib$N),"and",max(evsi$attrib$N),"as these are the boundaries within which the EVSI has been
                                                                                                           calculated. If the optimal sample size is calculated as either of these values it may be preferable to
                                                                                                           extend the limits of the EVSI calculation to find the optimal sample size."),
                                                                                                         p("Finally, note that the
                                                                                                           ENBS is a relatively \"flat\" function, meaning that, while the optimal sample size does exist, it is
                                                                                                           likely that sample sizes close to the optimal size will also give a similar benefit.",
                                                                                                           tags$br(),tags$br())),
                                                                                                fluidRow(column(5,p(h4("Optimal Sample Size: "),textOutput(outputId="SS"))),
                                                                                                         column(7,p(h4("Expected Net Benefit of Sampling: "),textOutput(outputId="ENBS"))),width=8))
                                                                                      )
                                                                                      ))
                                                                                    )
                           ),width=12
                           )
                           )
  server<-function(input,output,session){
    #save as output$hist - this will put hist in the hist part
    #inputs$n - access the inputs from the sliders ect...
    #Use inputs inside the render functions
    #inputs and outputs should be lists

    #render* functions work with the output functions to produce output
    #e.g. renderPlot({hist(rnorm(100))})

    #Dynamic Sliders
    output$PopDynam<-renderUI(sliderInput(
      inputId = "Pop",label="Incidence Population", max = input$Popmax, min = input$Popmin,value=c(input$Popmin,input$Popmax),step=1
    ))
    output$TimeDynam<-renderUI(sliderInput(
      inputId = "Time",label="Time Horizon", min = input$Timemin, max = input$Timemax,value=c(input$Timemin,input$Timemax),step=0.5
    ))
    output$Pop.OSDynam<-renderUI(sliderInput(
      inputId = "Pop.OS",label="Incidence Population", max = input$Popmax, min = input$Popmin,value=c(input$Popmax),step=1
    ))
    output$Time.OSDynam<-renderUI(sliderInput(
      inputId = "Time.OS",label="Time Horizon", min = input$Timemin, max = input$Timemax,value=c(input$Timemax),step=0.5
    ))

    #Reactive input min<max
    observeEvent(input$Timemin,{
      value.max<-max(input$Timemin,input$Timemax)
      updateNumericInput(session,"Timemax",value=value.max)
    })

    observeEvent(input$Timemax,{
      value.min<-min(input$Timemin,input$Timemax)
      updateNumericInput(session,"Timemin",value=value.min)
    })

    observeEvent(input$Popmin,{
      value.max<-max(input$Popmin,input$Popmax)
      updateNumericInput(session,"Popmax",value=value.max)
    })

    observeEvent(input$Popmax,{
      value.min<-min(input$Popmin,input$Popmax)
      updateNumericInput(session,"Popmin",value=value.min)
    })

    observeEvent(input$PerPersmin,{
      value.max<-max(input$PerPersmin,input$PerPersmax)
      updateNumericInput(session,"PerPersmax",value=value.max)
    })

    observeEvent(input$PerPersmax,{
      value.min<-min(input$PerPersmin,input$PerPersmax)
      updateNumericInput(session,"PerPersmin",value=value.min)
    })

    observeEvent(input$Setupmin,{
      value.max<-max(input$Setupmin,input$Setupmax)
      updateNumericInput(session,"Setupmax",value=value.max)
    })

    observeEvent(input$Setupmax,{
      value.min<-min(input$Setupmin,input$Setupmax)
      updateNumericInput(session,"Setupmin",value=value.min)
    })

    #Per person EVSI by Willingness to Pay - n slider
    output$ppEVSIbywtp<-renderPlot({
      N.chosen<-as.numeric(input$n)
      suppressWarnings(plot.evsi(evsi,N=N.chosen))
    }
    )

    #Per person EVSI by N - wtp slider
    output$ppEVSIbyn<-renderPlot({
      plot.evsi.N(evsi,wtp=as.numeric(input$wtp))
    })

    #Probability of Cost Effective Trial plot
    output$ProbCE<-renderPlot({
      if(is.null(input$Pop)){return(NULL)}
      if(is.null(input$Time)){return(NULL)}
      Pop<-as.numeric(input$Pop)
      Time<-as.numeric(input$Time)
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      N.chosen.CE<-evsi$attrib$N[which.min((evsi$attrib$N-as.numeric(input$n.CE))^2)]
      evsi.pop(evsi,setup=setup,pp=pp,
               Pop=Pop,Time=Time,Dis=input$Dis,
               wtp=as.numeric(input$wtp.CE),N=N.chosen.CE)

    })

    #Optimal Sample Size
    output$SS<-renderText({
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      suppressWarnings(optim.ss(evsi,setup,pp,input$Pop.OS,input$Time.OS,wtp=as.numeric(input$wtp.OS))$SS.max)
    }
    )

    output$ENBS<-renderText({
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      round(suppressWarnings(optim.ss(evsi,setup,pp,as.numeric(input$Pop.OS),as.numeric(input$Time.OS),
                                      wtp=as.numeric(input$wtp.OS))$ENBS),-1)
    })



  }
  shiny::shinyApp(ui=ui,server=server)
}
