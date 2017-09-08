##evsi################################################################
evsi<-function(model.stats,data,effects=NULL,costs=NULL,he=NULL,evi=NULL,parameters=NULL,Q=30,data.stats=NULL,update=c("bugs","jags"),
               n.burnin=1000,n.thin=1,n.iter=5000){
  ##'Compute the EVSI for a fixed design using the Heath et al. method
  ##'INPUTS:
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or JAGS model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. If the future data has already been generated then
  ##'   the data arguement can be given as a list. The elements of the data list should
  ##'   be data lists for JAGS or BUGS.
  ##'@param effects This can either be given as a string which defines the name of the
  ##'   effectivness measure in the BUGS/JAGS model. Or it can be given as a
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
  ##'@param data.stats A data file for the BUGS/JAGS model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the JAGS/BUGS model
  ##'@param n.burnin The burnin for the JAGS/BUGS model
  ##'@param n.iter The number of interations for the JAGS/BUGS model
  ##'
  ##OUTPUTS:
  ##'@return The EVSI for the specific design by willingness-to-pay. This can be plotted.
  ##'@examples NULL

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
  if(update=="jags"){ #Model can be written in jags.
    # Checks if rjags is installed (and if not, asks for it)
    if(!isTRUE(requireNamespace("rjags",quietly=TRUE))) {
      stop("You need to install the R package 'rjags' and the software 'JAGS'. \nPlease see http://mcmc-jags.sourceforge.net/ for instructions
           on installing 'JAGS' and then run in your R terminal:\n install.packages('rjags')")
    }
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
      #Runs the jags model based on the rjags package
      Model.JAGS<- rjags::jags.model(model.stats,data=data.stats,quiet=TRUE)
      update(Model.JAGS,n.burnin,progress.bar="none")
      Prior.Pred <- rjags::coda.samples(Model.JAGS, prior.pred.data, n.iter=n.iter,n.thin=n.thin,progress.bar="none")
      PP.sample<-as.data.frame(Prior.Pred[[1]])
    }

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
      # Checks if BCEA is installed (and if not, asks for it)
      if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }
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
      he<-BCEA::bcea(e,c)
    }
    #Check number of interventions
    #if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
    #                            multi-decision problems")}

    if(class(evi)!="evppi"){
      if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }
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
      evi<-BCEA::evppi(of.interest,PP.sample,he,h.value=0.00000001)
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

      Model.JAGS<- rjags::jags.model(model.stats,data =  data.full,quiet = TRUE)
      update(Model.JAGS,n.burnin,progress.bar="none")
      samples <- rjags::coda.samples(Model.JAGS, moniter, n.iter=n.iter,n.thin=n.thin,n.chain=1,progress.bar="none")
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
  if(update=="bugs"){#Model can be written in BUGS.
    # Checks if R2OpenBUGS is installed (and if not, asks for it)
    if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
      stop("You need to install the R package 'R2OpenBUGS' and the software 'OpenBUGS'. \nPlease see http://www.openbugs.net/w/FrontPage for instructions
           on installing 'OpenBUGS' and then run in your R terminal:\n install.packages('R2OpenBUGS')")
    }
    #Defines the model parameters that are required to calculate the costs and effects
    moniter<-names(which(sapply(nameofinterest,grep.fun,main=readLines(model.stats))>0))

    #Track all the variables of interest
    prior.pred.data<-unique(moniter)
    #Track both the data and the parameters of interest if either are needed.
    if(cl.dat=="character"){prior.pred.data <-unique(c(prior.pred.data,data))}
    if(length(parameters)!=0){prior.pred.data<-unique(c(prior.pred.data,parameters))}

    if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
      #Runs the BUGS model based on the R2OpenBUGS package
      Model.BUGS<- R2OpenBUGS::bugs(data.stats,inits=NULL,parameters.to.save=prior.pred.data,
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
      if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }
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
      he<-BCEA::bcea(e,c)
    }
    #Check number of interventions
    #if(he$n.comparisons>1){stop("WARNING:This EVSI calculation method is currently not implemented for
     #                           multi-decision problems")}

    if(class(evi)!="evppi"){
      if(!isTRUE(requireNamespace("BCEA",quietly=TRUE))) {
        stop("You need to install the R package 'BCEA'. Please run in your R terminal:\n install.packages('BCEA')")
      }
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
      evi<-BCEA::evppi(of.interest,PP.sample,he)
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

      Model.BUGS<- R2OpenBUGS::bugs(data.full,inits=NULL,parameters.to.save=moniter,
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

  CreateCov<-function(cov,i,j,d,k){
    cov.i.j<-k^2*cov[i,j]-k*(cov[i,d+j]+cov[d+i,j])+cov(d+i,d+j)
    return(cov.i.j)
  }

  #Calcualte the EVSI accross different WTP
  mean.var<-apply(simplify2array(var.prepost),1:2,mean)
  EVSI<-array(NA,dim=c(1,length(he$k),1))
  k.<-1
  #EVSI[,1:length(he$k),]<-foreach(k=he$k,.combine="c") %do% {
  for(k in he$k){
    #Fitted INB for each k
    INB<-as.matrix(k*evi$fitted.effects[,-he$comparators]-evi$fitted.costs[,-he$comparators])
    var.INB<-as.matrix(var(INB))
    var.full<-as.matrix(var(he$ib[k.,]))
    #Mean preposterior variance for each k
    Var<-matrix(NA, nrow=he$comparisons,ncol=he$comparisons)
    for(i in 1:he$comparisons){
      for(j in 1:i)
        Var[i,j]<-Var[j,i]<-CreateCov(mean.var,i,j,he$comparisons,k)
    }
    #Var<-k^2*mean.var[1,1]+mean.var[2,2]-2*k*mean.var[1,2]
    #Rescale fitted INB
    INB.star<-  sweep(as.matrix(sweep(INB,2,colMeans(INB)))%*%solve(expm::sqrtm(var.INB))%*%
                        expm::sqrtm(var.INB-Var),2,colMeans(INB),"+")
      #(INB-mean(INB))/sd(INB)*sqrt(max(0,var(he$ib[k.,])-Var))+mean(INB)
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

