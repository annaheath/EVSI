##comp.evsi.N############################################################
comp.evsi.N<-function(model.stats,data,N,N.range=c(30,1500),effects,costs,he=NULL,parameters=NULL,evi=NULL,Q=50,
                      data.stats=NULL,update=c("bugs","jags"),
                      n.burnin=1000,n.thin=1,n.iter=5000){
  ##'Function generates the data points that will be used to calculate the EVSI for different
  ##'sample sizes.
  ##INPUTS.
  ##'@param model.stats A .txt file containing the model file of a Bayesian model.
  ##'   This should be a BUGS or JAGS model.
  ##'@param data A string or vector of strings that defines the name of the future
  ##'   data in the model file. If the future data has already been generated then
  ##'   the data arguement can be given as a list. The elements of the data list should
  ##'   be data lists for JAGS or BUGS.
  ##'@param N A string defining the name of the variable controlling the sample size in the
  ##'   model.
  ##'@param N.range A two-vector defining the minimum and maximum values considered for the
  ##'   sample size. If the future data is given in the data arguement then N.range should
  ##'   give the sample sizes for which the data was generated.
  ##'@param effects This can either be given as a string which defines the name of the
  ##'   effectivness measure in the BUGS/JAGS model. Or it can be given as a
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
  ##'@param data.stats A data file for the BUGS/JAGS model. This is the data used to inform
  ##'   the base case analysis. If empty then it is assumed that the models are
  ##'@param update Defines the Bayesian engine that should be used to update the
  ##'   the Bayesian model file given in model.stats.
  ##'@param n.thin The thinning for the JAGS/BUGS model
  ##'@param n.burnin The burnin for the JAGS/BUGS model
  ##'@param n.iter The number of interations for the JAGS/BUGS model
  ##'
  ##'OUTPUTS.
  ##'@return An evsi.comp object to be used in the evsi.calc function.


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
    if(length(jags.params)!=0){
    multiple.params<-moniter[-jags.params]}
    else{multiple.params<-moniter}

    #Function to extract the matrix parameters from the jags object
    multiple.params.extract<-function(multi.params){
      if(length(multi.params)==0){return(NULL)}
      else{
        list.params<-list()
        length(list.params)<-length(multi.params)
        for(j in 1:length(multi.params)){
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

  quantiles<-sample((1:Q)/(Q+1),replace=F)
  if(length(N.range)==2){N.samp<-trunc((seq(sqrt(N.range[1]),sqrt(N.range[2]),length.out = Q))^2)}
  if(length(N.range)>2){N.samp<-N.range}
  while(abs(cor(quantiles,N.samp))>0.0005){quantiles<-sample(quantiles,Q,replace=FALSE)}
  
  #Generate future samples by finding prior-predictive distribution
  if(update=="jags"){#Model can be written in jags.
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
    if(class(evi)=="evppi"){parameters<-evi$index}
    prior.pred.data<-unique(c(prior.pred.data,parameters))


    ####Calculate EVSI by Quadrature####
    var.prepost<-list()
    start<-Sys.time()
    for(q in 1:Q){
      #Set the data including the sample size.
      Samp.Size<-list(N.samp[q])
      names(Samp.Size)<-N

      if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
        #Runs the jags model based on the rjags package
        Model.JAGS<- rjags::jags.model(model.stats,data=append(Samp.Size,data.stats),quiet=TRUE)
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
        index.params<-list()
        for(l in 1:length(parameters)){
          index.params[[l]]<-which(colnames(PP.sample)==parameters[l])
        }
        length.data<-length(unlist(index.data))
        names.data<-data
        Data.Fut<-array()
        #Selecting one row of the parameters of interest
        #Find the quantile of one of the parameters of interest, then select the row with that...
        index.pp<-sample(unlist(index.params),1)
        index.choose<-which(PP.sample[,index.pp]==quantile(PP.sample[,index.pp],probs=quantiles[q],type=3))
        for(d in 1:length.data){
          #Creating list of the future data to give to jags
          Data.Fut[unlist(index.data)[d]]<-as.numeric(PP.sample[index.choose,unlist(index.data)[d]])
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
        # Checks if BCEA is installed (and if not, asks for it)
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



      #Both data sets must be given to jags
      data.full<-append(append(Data.Fut,data.stats),Samp.Size)

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
      var.prepost[[q]]<-var(int)
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))
      }
      print(paste("Update",q,"completed"))
    }
    #Calcualte the EVSI accross different WTP
    print("Using curve fitting to find EVSI for alternative sample sizes")
    ##### GB: Need to define the variables var.PI and x, or else when compiling the package R will throw a message for no-bindings
    var.PI=x=NULL
    #####
    #Write the Model - remove need for R2OpenBUGS so can run with just JAGS.
    model.ab<-c("model
                {
                beta ~ dnorm(Nmax/2, shape.Nmax)  T(0.00000E+00, )
                for (i in 1:N) {
                y[i] ~ dnorm(mu[i], tau)
                mu[i] <- var.PI * (x[i]/(x[i] + beta))
                }
                sigma ~ dt(sigma.mu, sigma.tau, 3)  T(0.00000E+00, )
                tau <- 1/sigma^2
                }
                ")
    file.curve.fitting <- file.path("~",fileext="model_curve_fitting_EVSI.txt")
    writeLines(model.ab,file.curve.fitting)
    #######
    n.burnin <- 1000  # Number of burn in iterations
    n.thin<-1
    n.iter <- 3000 # Number of iterations per chain
    #######


    n.unique.entries<-he$n.comparisons*(he$n.comparisons+1)/2
    beta.mat<-array(NA,dim=c(3000,2,n.unique.entries))
    #for k over the wtp dimension
    e.var<-var(he$delta.e)
    c.var<-var(he$delta.c)
    e.fit<-evi$fitted.effects[,-he$n.comparators]
    e.fit.var<-var(e.fit)
    c.fit<-evi$fitted.costs[,-he$n.comparators]
    c.fit.var<-var(c.fit)

    n.entry<-1
    for(i in 1:he$n.comparisons){
      for(j in i:he$n.comparisons){
        if(he$n.comparisons>1){
          pre.var.e<-simplify2array(var.prepost)[1:he$n.comparisons,1:he$n.comparisons,]
          y<-t(e.var[i,j]-pre.var.e[i,j,])
        }
        if(he$n.comparisons==1){
          pre.var.e<-simplify2array(var.prepost)[1:he$n.comparisons,1:he$n.comparisons,]
          y<-e.var-pre.var.e
        }

        if(sd(y)==0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=100,
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }
        if(sd(y)!=0){
        data.a.b<-list(sigma.mu=sd(y)/2,
                       sigma.tau=1/sd(y),
                       N=length(N.samp),
                       shape.Nmax=0.0005/max(N.samp),
                       var.PI=as.matrix(e.fit.var)[i,j],
                       Nmax=max(N.samp),
                       y=as.vector(y),
                       x=as.vector(N.samp)
        )
        }
        Model.JAGS<- rjags::jags.model(file.curve.fitting,data=data.a.b,quiet=TRUE)
        update(Model.JAGS,n.burnin,progress.bar="none")
        beta.ab <- rjags::coda.samples(Model.JAGS, c("beta"), n.iter=n.iter,n.thin=n.thin,progress.bar="none")

        beta.mat[,1,n.entry]<-as.data.frame(beta.ab[[1]])[,1]
        n.entry<-n.entry+1
      }
    }

    n.entry<-1
    for(i in 1:he$n.comparisons){
      for(j in i:he$n.comparisons){
        if(he$n.comparisons>1){
          pre.var.c<-simplify2array(var.prepost)[(he$n.comparisons+1):(2*he$n.comparisons),
                                                 (he$n.comparisons+1):(2*he$n.comparisons),]
          y<-t(c.var[i,j]-pre.var.c[i,j,])
        }
        if(he$n.comparisons==1){
          pre.var.c<-simplify2array(var.prepost)[(he$n.comparisons+1):(2*he$n.comparisons),
                                                 (he$n.comparisons+1):(2*he$n.comparisons),]
          y<-c.var-pre.var.c
        }

        
        if(sd(y)==0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=100,
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }
        if(sd(y)!=0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=1/sd(y),
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }

        Model.JAGS<- rjags::jags.model(file.curve.fitting,data=data.a.b,quiet=TRUE)
        update(Model.JAGS,n.burnin,progress.bar="none")
        beta.ab <- rjags::coda.samples(Model.JAGS, c("beta"), n.iter=n.iter,n.thin=n.thin,progress.bar="none")

        beta.mat[,2,n.entry]<-as.data.frame(beta.ab[[1]])[,1]
        n.entry<-n.entry+1
      }
    }

  }

  #Analysis in BUGS
  if(update=="bugs"){ #Model can be written in BUGS.
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
    if(class(evi)=="evppi"){parameters<-evi$index}
    prior.pred.data<-unique(c(prior.pred.data,parameters))

    var.prepost<-list()
    start<-Sys.time()
    for(q in 1:Q){
      #Set the data including the sample size.
      Samp.Size<-list(N.samp[q])
      names(Samp.Size)<-N

      if((cl.dat=="character")|(class(he)!="bcea")|(class(evi)!="evppi")){
        #Runs the BUGS model based on the R2OpenBUGS package
        Model.BUGS<- R2OpenBUGS::bugs(append(data.stats,Samp.Size),inits=NULL,parameters.to.save=prior.pred.data,
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
        index.params<-list()
        for(l in 1:length(parameters)){
          index.params[[l]]<-which(colnames(PP.sample)==parameters[l])
        }
        length.data<-length(unlist(index.data))
        names.data<-data
        Data.Fut<-array()
        #Selecting one row of the parameters of interest
        #Find the quantile of one of the parameters of interest, then select the row with that...
        index.pp<-sample(unlist(index.params),1)
        index.choose<-which(PP.sample[,index.pp]==quantile(PP.sample[,index.pp],probs=quantiles[q],type=3))
        for(d in 1:length.data){
          #Creating list of the future data to give to jags
          Data.Fut[unlist(index.data)[d]]<-as.numeric(PP.sample[index.choose,unlist(index.data)[d]])
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
        evi<-BCEA::evppi(of.interest,PP.sample,he)
      }

      #Both data sets must be given to jags
      data.full<-append(append(Data.Fut,data.stats),Samp.Size)

      ####Calculate EVSI by Quadrature####
            Model.BUGS<- R2OpenBUGS::bugs(data.full,inits=NULL,parameters.to.save=moniter,
                        model.file=model.stats, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chains=1,
                        DIC=FALSE,debug=FALSE)
      #Create a dataframe containing all the BUGS samples
      sample<-as.data.frame(Model.BUGS$sims.matrix)

      #Use the BUGS output as inputs for the effects and costs functions
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
      var.prepost[[q]]<-var(int)
      if(q==1){end<-Sys.time()
      comp.time<-end-start
      print(paste(c("Model updating requires ",round(comp.time,2)," seconds. The EVSI will be calculated using ",Q," model updates. The remaining computation time is around ",round(comp.time*(Q-1)/60,0)," minutes. The current time is ",strftime(Sys.time())),
                  sep="",collapse = ""))
      }
      print(paste("Update",q,"completed"))
    }
    #Calcualte the EVSI accross different WTP
    print("Using curve fitting to find EVSI for alternative sample sizes")
    ##### GB: Need to define the variables var.PI and x, or else when compiling the package R will throw a message for no-bindings
    var.PI=x=NULL
    #####
    #Write the Model
    model.ab<-c("model
                {
                beta ~ dnorm(Nmax, shape.Nmax)  I(0.00000E+00, )
                for (i in 1:N) {
                y[i] ~ dnorm(mu[i], tau)
                mu[i] <- var.PI * (x[i]/(x[i] + beta))
                }
                sigma ~ dt(sigma.mu, sigma.tau, 3)  I(0.00000E+00, )
                tau <- 1/sigma*sigma
                }
                ")
    file.curve.fitting <- file.path(getwd(),fileext="model_curve_fitting_EVSI.txt")
    writeLines(model.ab,file.curve.fitting)

    #######
    n.burnin <- 1000  # Number of burn in iterations
    n.thin<-1
    n.iter <- 3000 # Number of iterations per chain
    #######


    n.unique.entries<-he$n.comparisons*(he$n.comparisons+1)/2
    beta.mat<-array(NA,dim=c(3000,2,n.unique.entries))
    #for k over the wtp dimension
    e.var<-var(he$delta.e)
    c.var<-var(he$delta.c)
    e.fit<-evi$fitted.effects[,-he$n.comparators]
    e.fit.var<-var(e.fit)
    c.fit<-evi$fitted.costs[,-he$n.comparators]
    c.fit.var<-var(c.fit)

    for(i in 1:he$n.comparisons){
      for(j in i:he$n.comparisons){
        if(he$n.comparisons>1){
        pre.var.e<-simplify2array(var.prepost)[1:he$n.comparisons,1:he$n.comparisons,]
        y<-t(e.var[i,j]-pre.var.e[i,j,])
        }
        if(he$n.comparisons==1){
          pre.var.e<-simplify2array(var.prepost)[1:he$n.comparisons,1:he$n.comparisons,]
          y<-e.var-pre.var.e
        }

        y<-t(e.var[i,j]-pre.var.e[i,j])

        
        if(sd(y)==0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=100,
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }
        if(sd(y)!=0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=1/sd(y),
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }

        beta.ab<- R2OpenBUGS::bugs(data.a.b,inits=NULL,parameters.to.save=c("beta"),
                                      model.file=file.curve.fitting, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chains=1,
                                      DIC=FALSE,debug=FALSE)

        beta.mat[,1,n.entry]<-as.data.frame(beta.ab$sims.matrix)[,1]
        n.entry<-n.entry+1
      }
    }

    n.entry<-1
    for(i in 1:he$n.comparisons){
      for(j in i:he$n.comparisons){
        if(he$n.comparisons>1){
          pre.var.c<-simplify2array(var.prepost)[(he$n.comparisons+1):(2*he$n.comparisons),
                                                 (he$n.comparisons+1):(2*he$n.comparisons),]
          y<-t(c.var[i,j]-pre.var.c[i,j,])
        }
        if(he$n.comparisons==1){
          pre.var.c<-simplify2array(var.prepost)[(he$n.comparisons+1):(2*he$n.comparisons),
                                                 (he$n.comparisons+1):(2*he$n.comparisons),]
          y<-c.var-pre.var.c
        }

        
        if(sd(y)==0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=100,
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }
        if(sd(y)!=0){
          data.a.b<-list(sigma.mu=sd(y)/2,
                         sigma.tau=1/sd(y),
                         N=length(N.samp),
                         shape.Nmax=0.0005/max(N.samp),
                         var.PI=as.matrix(e.fit.var)[i,j],
                         Nmax=max(N.samp),
                         y=as.vector(y),
                         x=as.vector(N.samp)
          )
        }
        beta.ab<- R2OpenBUGS::bugs(data.a.b,inits=NULL,parameters.to.save=c("beta"),
                                   model.file=file.curve.fitting, n.burnin=n.burnin,n.iter=n.iter+n.burnin,n.thin=n.thin,n.chains=1,
                                   DIC=FALSE,debug=FALSE)
        
        beta.mat[,1,n.entry]<-as.data.frame(beta.ab$sims.matrix)[,1]
        n.entry<-n.entry+1
      }
    }
  }


  #Return EVSI, plus evppi object and bcea object to plot EVSI.
  to.return<-list(beta = beta.mat,
                  N=N.samp,
                  evppi=evi,
                  he=he)
  class(to.return)<-"evsi.N"
  return(to.return)
}
