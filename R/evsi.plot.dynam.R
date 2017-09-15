##evsi.plot.dynam######################################################
evsi.plot.dynam<-function(evsi,width=600,height=600){
  ##'Function launches shiny app to show dynamic graphics
  ##'INPUTS
  ##'@param evsi An evsi object for which dynamic graphics are needed.
  ##'@param width The width of the graph produced in the webapp (default 600pt)
  ##'@param height The height of the graph produced in the webapp (default 600pt)
  ##'
  ##OUTPUTS
  ##'@return A dynamic shiny app that plots the EVSI for different wtp,
  ##'sample size, incidence population, time horizon and costs as required.

  # Checks if shiny and shinythemes are installed (and if not, asks for them)
  if(!isTRUE(requireNamespace("shiny",quietly=TRUE))) {
    stop("You need to install the R package 'shiny'.Please run in your R terminal:\n install.packages('shiny')")
  }
  if(!isTRUE(requireNamespace("shinythemes",quietly=TRUE))) {
    stop("You need to install the R package 'shiny'.Please run in your R terminal:\n install.packages('shinythemes')")
  }

  digits<-nchar(as.character(round(diff(evsi$attrib$wtp))))[1]
  if(is.na(digits)){rounding<-0}
  else{rounding<-1-digits}

  #if(class(evsi$attrib$N)!="numeric"){
  #  evsi$attrib$N<-0
  #}
  ui<-shiny::fluidPage(theme=shinythemes::shinytheme('united'),
                       shiny::titlePanel(shiny::h1("Visualisations for the EVSI")),
                       shiny::mainPanel(
                         #Create Tabs
                         shiny::tabsetPanel(
                           #Description Panel
                           shiny::tabPanel("Introduction to plotting",
                                           shiny::fluidRow(shiny::p("The Expected Value of Sample Information can be used to determine the economic benefit of a future
                                               trial. However, in general, this economic benefit is conditional on several deterministic inputs that can
                                               generally have a large impact on the value of the EVSI and any decision you make with the EVSI in mind.
                                               Therefore, we have developed this tool to allow for the simple visualistion of the EVSI along with easy
                                               manipulation of these deterministic inputs. This tool also allows for a simple user interface that can be
                                               shown to collaborators and other stakeholders.",
                                                                    shiny::h2("Willingness to Pay"),
                                                                    shiny::p("The willingness to pay (WTP) is the amount of money the decision maker has avaliable to pay for an
                                                 additional 1 unit of benefit. This is typically given as a range and therefore the EVSI can be visualised
                                                 across different WTP values in the tab WTP. For the other analyses the WTP must be fixed. Although the plots
                                                 can easily be redrawn for different values of the WTP allowing the analysis to be easily completed for a
                                                 large number of different thresholds."),
                                                                    shiny::h2("Sample Size"),
                                                                    shiny::p("The sample size of the future trial is rarely determined before the EVSI analysis and one of the analysis
                                                 avaliable in this tool involves determining the optimal sample size for the future trial. The EVSI can be
                                                 visualised by sample size on the tab N. The value of a sample is bounded above by the value of resolving
                                                 all uncertainty in the parameters of interest. The key information to be gleaned from the by N plot is the
                                                 speed at which the EVSI reachs this upper bound."),
                                                                    shiny::h2("Trial Cost-effectivenes"),
                                                                    shiny::p("A key use of the EVSI is to determine whether a future trial will be cost-effective. This means that the
                                                 value of the future trial exceeds the cost of the trial. While in general this is a simple extension of the
                                                 ideas underpinning the EVSI, it typically depends on some additional inputs which are rarely known with
                                                 certainty."))),
                                           shiny::fluidRow(shiny::column(11,offset=1,shiny::h3("Trial costs"),
                                                                         shiny::p("The cost of the trial must be specified to allow comparison with the EVSI. However, these are rarely known
                                                      with certainty and so this application allows a range of values to be specified for the cost and then
                                                      considers the cost-effectiveness of the trial taking into account this uncertainty."),
                                                                         shiny::h3("Incidence Population"),
                                                                         shiny::p("The EVSI is calculated as the EVSI per person who will benefit from the treatments under consideration.
                                                      Therefore, to compare with the trial costs, we must multiply by the number of people who will benefit
                                                      from the treatments. We call this the \"Incidence Population\". While there may be some literature that
                                                      can inform this parameter it cannot be known and therefore the cost-effectiveness of the trial is determined
                                                      for different levels of this population. This is an important consideration as the level of this population
                                                      can make a large difference to the cost-effectiveness of the trial."),
                                                                         shiny::h3("Time Horizon"),
                                                                         shiny::p("Finally, the time horizon of a treatment is the length of time the treatments will be available in the
                                                      market before an alternative superseeds the treatments by being more cost-effective. Typically, this is
                                                      assumed to be around 10 years but we allow variation in this parameter to consider the cost-effectiveness
                                                      for alternative time horizons. This is powerful as it is rarely known the length of time a technology will
                                                      be available in the market.")))
                                                    ),

                           #Per person EVSI by Willingness to Pay - n slider
                           shiny::tabPanel("EVSI by Willingness To Pay",
                                           shiny::sidebarPanel(shiny::p("The EVSI changes depending on the willingness-to-pay of the underlying decision maker and therefore
                                                   it is useful to visualise the EVSI for different WTP thresholds. This plot is also used visualise the
                                                   EVPI - which gives the maximum value for ANY future trial - and the EVPPI - which gives a maximum
                                                   value for a study targeting uncertainty in the parameters of interest of our study.
                                                   The EVSI also changes for different sample sizes as studies increase in value as the sample size increases.
                                                   The sample size can therefore be changed:"),
                                                               shiny::selectInput(inputId="n",label="Choose a sample size",
                                                             choices=evsi$attrib$N,
                                                             selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)]),
                                                             shiny::p("The EVSI should remain below the EVPPI for all sample sizes and willingness-to-pay thresholds.
                                                   All three measures should reach a sharp peak - this represents the \"break-even\" point between the
                                                   two treatment options where the decision uncertainty is at its maximum as the two treatment options are
                                                   equally likely to be cost-effective."),width=4),
                                           shiny::mainPanel(shiny::plotOutput(outputId="ppEVSIbywtp",width=width,height=height),
                                                            shiny::p("Citations: McCabe, C., Claxton K. and Culyer, A., ",
                                                                     shiny::a("The NICE cost-effectiveness threshold",href="https://link.springer.com/article/10.2165/00019053-200826090-00004"),
                                                                     shiny::em(", PharmacoEconomics"),", 2008"),
                                                            shiny::p("Baio, G., ",
                                                                     shiny::a("Bayesian Methods in Health Economics",href="https://sites.google.com/a/statistica.it/gianluca/bookhe"),
                                                                     shiny::em(", Springer"),", 2012"),
                                                            shiny::p("Baio, G., Berardi, A. and Heath, A., ",
                                                                     shiny::a("Bayesian Cost-Effectiveness Analysis with the R package BCEA",href="http://www.springer.com/gb/book/9783319557168"),
                                                                     shiny::em(", Springer"),", 2017"),width=8)
                                                 ),

                           #Per person EVSI by N - wtp slider

                           shiny::tabPanel("EVSI by Sample Size",
                                           shiny::sidebarPanel(shiny::p("The EVSI increases as the sample size of the underlying trial increases. This graphic shows the
                                                   EVSI across different sample sizes. This relationship with N changes depending on the value of the
                                                   willingness-to-pay so the plot can be considered for changing values of the WTP.
                                                   The EVSI is calculated using Bayesian regression. This means that posterior credible intervals for the EVSI
                                                   can be calculated and these are plotted on the graphic. This demonstrates the uncertainty in the EVSI
                                                   estimate. If the uncertainty is too large for decision making then Q should be increased in the
                                                   EVSI calculation."),
                                                               shiny::selectInput(inputId="wtp",label="Choose a Willingness-to-Pay Threshold",
                                                             choices=round(evsi$attrib$wtp,rounding),
                                                             selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2)],rounding)),
                                                 #sliderInput(inputId="wtp",label="Choose a Willingness-to-Pay Threshold",
                                                 #min=min(evsi$attrib$wtp),max=max(evsi$attrib$wtp),value=min(evsi$he$kstar),
                                                 #step=round(evsi$attrib$wtp[2]-evsi$attrib$wtp[1])),
                                                 shiny::p(paste("In general, the EVSI will be more accurately estimated for higher values of the EVSI. This
                                                         would will be for willingness to pay values close to the \"break-even\" point of ",
                                                         round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2)],rounding),".")),width=4),
                                           shiny::mainPanel(shiny::plotOutput(outputId="ppEVSIbyn",width=width,height=height),
                                                            shiny::p("Citations: Heath, A., Manolopoulou I. and Baio, G., ",
                                                shiny::a("Efficient Monte Carlo Estimation of the Expected Value of Sample Information using Moment Matching",href="https://arxiv.org/abs/1611.01373"),
                                                shiny::em(", arXiv Preprint"),", 2017"),
                                                shiny::p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                         shiny::a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                         shiny::em(", Working Paper"),", 2017"),
                                              width=8)
                                    ),

                           #Trial Costs Input
                           shiny::tabPanel("Cost-effectiveness of a Trial",shiny::column(12,
                                                                                         shiny::tabsetPanel(id="CE",shiny::tabPanel("Trial Costs",shiny::fluidRow(
                                                                                           shiny::column(5,
                                                                                                         shiny::fluidRow(shiny::column(12,shiny::p("To determine the cost-effectiveness of a trial, the EVSI must be compared with the costs of
                                                                                                         undertaking the trial. In general, these costs are split into two categories; the setup and the per person costs.
                                                                                                         The are setup costs are the overhead costs of the trial and will be incurred irrespective of the size
                                                                                                         of the trial. Typical setup costs may be training for staff or the purchase of specialised equipment.
                                                                                                         Additional costs will then be incurred by each participant in the trial such as the cost of
                                                                                                         administering the treatment or following up the patient."),
                                                                                                                                       shiny::p("It is unlikely that these costs will be known exactly. Therefore, you should give a range of possible
                                                                                                      values for the costs. If the costs are known with certainty then simply input the same values for
                                                                                                      the maximum and minimum possible values of the costs. ")
                                                                                                    )),
                                                                                                    shiny::fluidRow(
                                                                                                      shiny::column(6,shiny::numericInput(inputId="Setupmin",label="Minimum Setup Costs for the Trial",
                                                                                                            value=1000,step=10,min=0),
                                                                                                            shiny::numericInput(inputId="PerPersmin",label="Minimum Cost Per Person",min=0,
                                                                                                          value=1500,step=10)),
                                                                                                      shiny::column(6,shiny::numericInput(inputId="Setupmax",label="Maximum Setup Costs for the Trial",
                                                                                                            value=100,step=10,min=0),
                                                                                                            shiny::numericInput(inputId="PerPersmax",label="Maximum Cost Per Person",min=0,
                                                                                                          value=1000,step=10)
                                                                                      ))),
                                                                                      shiny::column(5,shiny::fluidRow(shiny::p("The cost-effectiveness of a trial depends on two additional inputs. These are the incidence population, i.e.
                                                                                                 the yearly incidence of the disease under consideration. This gives the number of patients that will benefit from
                                                                                                 the treatment in each year that the treatment is used. This may not be known with certainty and the plot considers
                                                                                                 the cost-effectiveness of the trial for different possible numbers of patients. However, it is necessary to give
                                                                                                 possible values of the incidence population."),
                                                                                                                      shiny::p("The time horizon gives the number of years that the most cost-effective treatment will be available.
                                                                                                 This can be thought of as the number of years before a more effective treatment will be developed.
                                                                                                 This will depend on the disease areas as fast moving diseases such as cancer will have a shorter time
                                                                                                 horizon. The maximum and minimum possible values for the time horizon should be specified here.")),
                                                                                                    shiny::fluidRow(shiny::column(6,shiny::numericInput(inputId="Popmin",label="Minimum Incidence Population",value=0,step=100,min=0),
                                                                                                                                  shiny::numericInput(inputId="Timemin",label="Minimum Time Horizon",value=0,step=1,min=0)),
                                                                                                                    shiny::column(6,
                                                                                                                                  shiny::numericInput(inputId="Popmax",label="Maximum Incidence Population",value=1e+05,step=100,min=0),
                                                                                                                                  shiny::numericInput(inputId="Timemax",label="Maximum Time Horizon",min=0,value=25,step=1)
                                                                                             ))),
                                                                                      shiny::column(2,shiny::p("The final input to determine the cost-effectiveness of the trial is the discount rate for future treatments.
                                                                                        In general, health benefits now are more valuable than health benefits in the future. NICE recommend 3.5% as the
                                                                                        discount rate for the treatments but this can be changed here."),
                                                                                                    shiny::numericInput(inputId="Dis",label="Discount Rate",
                                                                                                 value=0.035,step=0.001))),
                                                                                      shiny::fluidRow(shiny::column(12,shiny::p("Citations:")))),
                                                                                      shiny::tabPanel("Probability of CE Trial",
                                                                                                      shiny::fluidRow(shiny::sidebarPanel(#Population EVSI - prob of CE plot
                                                                                                        shiny::selectInput(inputId="n.CE",label="Choose a sample size",
                                                                                                    choices=evsi$attrib$N,
                                                                                                    selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)]),
                                                                                                    shiny::selectInput(inputId="wtp.CE",label="Choose a Willingness-to-Pay Threshold",
                                                                                                    choices=round(evsi$attrib$wtp,rounding),selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2)],rounding)),
                                                                                                    shiny::uiOutput("PopDynam"),
                                                                                                    shiny::uiOutput("TimeDynam"),width=4),
                                                                                                    shiny::mainPanel(shiny::plotOutput(outputId="ProbCE",width=width,height=height),
                                                                                                  width=8)),
                                                                                                  shiny::fluidRow(shiny::column(12,shiny::p("The probability of having a cost effective trial is shown for different time horizons and incidence
                                                                                                           populations. The plot is",shiny::em("white"),"when the trial is cost-effective and",shiny::em("black"),"when the trial
                                                                                                           is not cost-effective. A trial is more likely to be cost effective as the time horizon and incidence
                                                                                                           population increases and therefore the plot will typically be white in the top right hand corner."),
                                                                                                                                shiny::p("A trial is",shiny::em("cost-effective"),"when the expected value of the information gained in the trial
                                                                                                        is greater than the trial costs. As both the trial costs and the EVSI are not known with certainty
                                                                                                        there are some values for the time horizon and incidence population where the costs of the trial could
                                                                                                        exceed the value but it is uncertain whether this will be the case. In these settings, the plot is
                                                                                                        blue. The darker the blue, the more likely it is that the costs of the trial exceed the EVSI."),
                                                                                                                                shiny::p("Citations: Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                                                         shiny::a("EVSI Visualisations for cost-effective trial analysis",href="https://arxiv.org/abs/1611.01373"),
                                                                                                                                         shiny::em(", NOT WRITTEN"),", 2018"),
                                                                                                                                shiny::p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                                                         shiny::a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                                                                                                         shiny::em(", Working Paper"),", 2017")))
                                                                                      ),
                                                                             #Optimal Sample Size

                                                                             shiny::tabPanel("Optimal Sample Size",value = "OSS",
                                                                                             shiny::sidebarPanel(shiny::selectInput(inputId="wtp.OS",label="Choose a Willingness-to-Pay Threshold",
                                                                                                               choices=round(evsi$attrib$wtp,rounding),
                                                                                                               selected=round(evsi$attrib$wtp[which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2)],rounding)),
                                                                                                               shiny::uiOutput("Pop.OSDynam"),
                                                                                                               shiny::uiOutput("Time.OSDynam"),
                                                                                                               shiny::p("It is possible to find the sample size for your trial that will give the maximum value for money.
                                                                                                                        In some cases this optimal sample size will be less than 0. In these cases, the EVSI indicates that the trial
                                                                                                                        is not cost-effective and the ",shiny::em("Expected Net Benefit of Sampling"),"(ENBS) is less than 0."),
                                                                                                               shiny::p("Note that the optimal sample size can only be found between",min(evsi$attrib$N),"and",max(evsi$attrib$N),
                                                                                                                        "as these are the boundaries within which the EVSI has been calculated. If the optimal sample size is given
                                                                                                                        as either of these values you will need to recalculate the EVSI for alternative values of N to find the true
                                                                                                                        optimal sample size."),
                                                                                                               shiny::p("Finally, the plot of the ENBS may demonstrate that there are large number of sample sizes for which
                                                                                                                        the ENBS is close to the optimal - any of these samples will give a similar level of benefit.")
                                                                                                               ,width=4),

                                                                                             shiny::mainPanel(shiny::fluidRow(shiny::column(5,shiny::p(shiny::h4("Optimal Sample Size: "),shiny::textOutput(outputId="SS"))),
                                                                                                                           shiny::column(7,shiny::p(shiny::h4("Expected Net Benefit of Sampling: "),shiny::textOutput(outputId="ENBS"))),
                                                                                                                           shiny::plotOutput("ENBS.plot",width=width,height=height),width=8))
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
    output$PopDynam<-shiny::renderUI(shiny::sliderInput(
      inputId = "Pop",label="Incidence Population", max = input$Popmax, min = input$Popmin,value=c(input$Popmin,input$Popmax),step=1
    ))
    output$TimeDynam<-shiny::renderUI(shiny::sliderInput(
      inputId = "Time",label="Time Horizon", min = input$Timemin, max = input$Timemax,value=c(input$Timemin,input$Timemax),step=0.5
    ))
    output$Pop.OSDynam<-shiny::renderUI(shiny::sliderInput(
      inputId = "Pop.OS",label="Incidence Population", max = input$Popmax, min = input$Popmin,value=c(input$Popmax),step=1
    ))
    output$Time.OSDynam<-shiny::renderUI(shiny::sliderInput(
      inputId = "Time.OS",label="Time Horizon", min = input$Timemin, max = input$Timemax,value=c(input$Timemax),step=0.5
    ))

    #Reactive input min<max
    shiny::observeEvent(input$Timemin,{
      value.max<-max(input$Timemin,input$Timemax)
      shiny::updateNumericInput(session,"Timemax",value=value.max)
    })

    shiny::observeEvent(input$Timemax,{
      value.min<-min(input$Timemin,input$Timemax)
      shiny::updateNumericInput(session,"Timemin",value=value.min)
    })

    shiny::observeEvent(input$Popmin,{
      value.max<-max(input$Popmin,input$Popmax)
      shiny::updateNumericInput(session,"Popmax",value=value.max)
    })

    shiny::observeEvent(input$Popmax,{
      value.min<-min(input$Popmin,input$Popmax)
      shiny::updateNumericInput(session,"Popmin",value=value.min)
    })

    shiny::observeEvent(input$PerPersmin,{
      value.max<-max(input$PerPersmin,input$PerPersmax)
      shiny::updateNumericInput(session,"PerPersmax",value=value.max)
    })

    shiny::observeEvent(input$PerPersmax,{
      value.min<-min(input$PerPersmin,input$PerPersmax)
      shiny::updateNumericInput(session,"PerPersmin",value=value.min)
    })

    shiny::observeEvent(input$Setupmin,{
      value.max<-max(input$Setupmin,input$Setupmax)
      shiny::updateNumericInput(session,"Setupmax",value=value.max)
    })

    shiny::observeEvent(input$Setupmax,{
      value.min<-min(input$Setupmin,input$Setupmax)
      shiny::updateNumericInput(session,"Setupmin",value=value.min)
    })

    #Per person EVSI by Willingness to Pay - n slider
    output$ppEVSIbywtp<-shiny::renderPlot({
      if(length(evsi$attrib$N)==1){
        suppressWarnings(plot.evsi(evsi))
      }
      if(length(evsi$attrib$N)>1){
      N.chosen<-as.numeric(input$n)
      suppressWarnings(plot.evsi(evsi,N=N.chosen))
}
    }
    )

    #Per person EVSI by N - wtp slider
    output$ppEVSIbyn<-shiny::renderPlot({
      plot.evsi.N(evsi,wtp=as.numeric(input$wtp))
    })

    #Probability of Cost Effective Trial plot
    output$ProbCE<-shiny::renderPlot({
      if(is.null(input$Pop)){return(NULL)}
      if(is.null(input$Time)){return(NULL)}
      Pop<-as.numeric(input$Pop)
      Time<-as.numeric(input$Time)
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      if(length(evsi$attrib$N)==1){
        evsi.pop(evsi,setup=setup,pp=pp,
                 Pop=Pop,Time=Time,Dis=input$Dis,
                 wtp=as.numeric(input$wtp.CE))
      }
      if(length(evsi$attrib$N)>1){
        N.chosen.CE<-evsi$attrib$N[which.min((evsi$attrib$N-as.numeric(input$n.CE))^2)]
        evsi.pop(evsi,setup=setup,pp=pp,
                 Pop=Pop,Time=Time,Dis=input$Dis,
                 wtp=as.numeric(input$wtp.CE),N=N.chosen.CE)
      }




    })

    #Optimal Sample Size
    output$SS<-shiny::renderText({
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      if(length(evsi$attrib$N)==1){return("Only one sample size considered")}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      suppressWarnings(optim.ss(evsi,setup,pp,input$Pop.OS,input$Time.OS,Dis=input$Dis,wtp=as.numeric(input$wtp.OS))$SS.max)
    }
    )

    output$ENBS<-shiny::renderText({
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      if(length(evsi$attrib$N)==1){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      round(suppressWarnings(optim.ss(evsi,setup,pp,as.numeric(input$Pop.OS),as.numeric(input$Time.OS),Dis=input$Dis,
                                      wtp=as.numeric(input$wtp.OS))$ENBS),-1)
    })

    output$ENBS.plot<-shiny::renderPlot({
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      if(length(evsi$attrib$N)==1){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))

      enbs.plot(evsi,setup,pp,Pop=as.numeric(input$Pop.OS),Time=as.numeric(input$Time.OS),
                Dis=input$Dis,wtp=as.numeric(input$wtp.OS))
    })



  }
  shiny::shinyApp(ui=ui,server=server)
}
