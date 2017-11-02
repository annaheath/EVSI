##launch.App######################################################
launch.App<-function(...){
  ##'Function launches shiny app to show dynamic graphics
  ##'INPUTS
  ##'@param ... Either empty or an evsi object. If given an evsi object then graphics will be displayed
  ##'           If empty then you can upload the EVSI using an upload page.
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

  if(length(list(...))>0){
    if(class(...)=="evsi"){
      simple<-"y"
      complex<-"n"
      obj<-list(...)
      evsi<-function(evsi=obj){
        evsi<-evsi[[1]]
        return(evsi)
      }
    }
    if(class(...)!="evsi"){
      simple<-"n"
      complex<-"y"
    }
  }

  if(length(list(...))==0){
    simple<-"n"
    complex<-"y"
  }

  ui<-shiny::fluidPage(theme=shinythemes::shinytheme('united'),
                       shiny::titlePanel("Visualisations for the EVSI"),
                       #Preserve proportionality of graphics
                       shiny::tags$head(shiny::tags$style(".shiny-plot-output{height:65vh !important;}")),
                       shiny::tags$head(shiny::tags$style(".shiny-plot-output{width:65vh !important;}")),
                       shiny::mainPanel(shiny::uiOutput("main"),width=12))

  server<-function(input,output,session){
    shiny::observe({
      inFile2<-input$inputs
      if(is.null(inFile2)) {return(NULL)}
      param = read.csv(inFile2$datapath, sep=',')
      updateSelectInput(session, 'parameter',choices = names(param))
    })

    inputs<-shiny::reactive({
      inFile2<-input$inputs
      if(is.null(inFile2)) {return(NULL)}
      inputs = read.csv(inFile2$datapath, sep=',')
    })

    effs.costs<-shiny::reactive({
      inFile1<-input$ec
      if(is.null(inFile1)){return(NULL)}
      ec.mat<-read.csv(inFile1$datapath, sep=',')
    })

    EVSI.mat<-shiny::reactive({
      inFile3<-input$evsi.csv
      if(is.null(inFile3)){return(NULL)}
      EVSI.mat<-read.csv(inFile3$datapath, sep=',',header=FALSE)
    })

    if(complex=="y"){evsi<-shiny::reactive({

      shiny::req(input$RunModel)
      shiny::req(effs.costs())
      shiny::req(EVSI.mat())
      shiny::isolate({
        N<-as.numeric(unlist(strsplit(input$N.input,",")))
        if(length(N)==0){
          N<-NULL
        }
        wtp<-as.numeric(unlist(strsplit(input$wtp.input,",")))

        if(length(wtp)==0){
          wtp<-NULL
        }
        parameter<-input$parameter

        n.cols<-as.integer(dim(effs.costs())[[2]])
        odds<-seq(1,n.cols,by=2)
        evens<-seq(2,n.cols,by=2)

        e<-as.matrix(effs.costs()[,odds])
        c<-as.matrix(effs.costs()[,evens])

        input.mat<-inputs()

        evsi.<-EVSI::evsi.upload(e,c,parameter,input.mat,EVSI.mat=as.matrix(EVSI.mat()),wtp=wtp,N=N)
      })
    })
    }


    shiny::observeEvent(input$RunModel,{
      output$BCEA<-shiny::renderPlot({
        evsi<-evsi()
        BCEA::plot.bcea(evsi$he)
      })
    })

    if(complex=="y"){shiny::observeEvent(input$Graphs,{
      evsi<-evsi()

      digits<-nchar(as.character(round(diff(evsi$attrib$wtp))))[1]
      if(is.na(digits)){rounding<-0}
      else{rounding<-1-digits}

      #Update inputs based on the EVSI upload
      shiny::updateSelectInput(session,"n",choices=evsi$attrib$N,selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)])
      shiny::updateSelectInput(session,"wtp",choices=round(evsi$attrib$wtp,rounding),
                               selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))
      shiny::updateSelectInput(session, "n.CE",choices=evsi$attrib$N,
                               selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)])
      shiny::updateSelectInput(session, "wtp.CE",choices=round(evsi$attrib$wtp,rounding),selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))
      shiny::updateSelectInput(session, "wtp.OS",choices=round(evsi$attrib$wtp,rounding),
                               selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))

    })
    }

    if(simple=="y"){

      evsi<-evsi()

      digits<-nchar(as.character(round(diff(evsi$attrib$wtp))))[1]
      if(is.na(digits)){rounding<-0}
      else{rounding<-1-digits}

      #Update inputs based on the EVSI upload
      shiny::updateSelectInput(session,"n",choices=evsi$attrib$N,selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)])
      shiny::updateSelectInput(session,"wtp",choices=round(evsi$attrib$wtp,rounding),
                               selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))
      shiny::updateSelectInput(session, "n.CE",choices=evsi$attrib$N,
                               selected=evsi$attrib$N[round(length(evsi$attrib$N)/2)])
      shiny::updateSelectInput(session, "wtp.CE",choices=round(evsi$attrib$wtp,rounding),selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))
      shiny::updateSelectInput(session, "wtp.OS",choices=round(evsi$attrib$wtp,rounding),
                               selected=round(evsi$attrib$wtp[max(1, which.min((evsi$attrib$wtp-evsi$he$kstar[1])^2))],rounding))


    }

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
      evsi<-evsi()
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
      evsi<-evsi()
      plot.evsi.N(evsi,wtp=as.numeric(input$wtp))
    })

    #Probability of Cost Effective Trial plot
    output$ProbCE<-shiny::renderPlot({
      evsi<-evsi()
      if(is.null(input$Pop)){return(NULL)}
      if(is.null(input$Time)){return(NULL)}
      Pop<-as.numeric(input$Pop)
      Time<-as.numeric(input$Time)
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      if(length(evsi$attrib$N)==1){
        evsi.pop(evsi,trial.cost=setup,
                 Pop=Pop,Time=Time,Dis=input$Dis,
                 wtp=as.numeric(input$wtp.CE))
      }
      if(length(evsi$attrib$N)>1){
        pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
        N.chosen.CE<-evsi$attrib$N[which.min((evsi$attrib$N-as.numeric(input$n.CE))^2)]
        evsi.pop(evsi,setup=setup,pp=pp,
                 Pop=Pop,Time=Time,Dis=input$Dis,
                 wtp=as.numeric(input$wtp.CE),N=N.chosen.CE)
      }




    })

    #Optimal Sample Size
    output$SS<-shiny::renderText({
      evsi<-evsi()
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))
      suppressWarnings(optimal<-optim.ss(evsi,setup,pp,input$Pop.OS,input$Time.OS,Dis=input$Dis,wtp=as.numeric(input$wtp.OS)))
      paste("The optimal sample size for this study is marked by a red triangle on the graph above and is equal to ",
            optimal$SS.max,
            ". However, any study with the sample size between",
            optimal$SS.I[1]," and ",optimal$SS.I[2],
            " has a value within 5% of this optimal value - where this area is marked with the red line on the graph.",sep="")
      
      
    }
    )

    #output$SS.min<-shiny::renderText({
    #  evsi<-evsi()
    #  if(is.null(input$Pop.OS)){return(NULL)}
    #  if(is.null(input$Time.OS)){return(NULL)}
    #  pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
    #  setup<-as.numeric(c(input$Setupmin,input$Setupmax))
    #  suppressWarnings(optimal<-optim.ss(evsi,setup,pp,input$Pop.OS,input$Time.OS,Dis=input$Dis,wtp=as.numeric(input$wtp.OS)))
    #  optimal$SS.max
    #  optimal$SS.I[1]
    #}
    #)
    #output$SS.max<-shiny::renderText({
    #  evsi<-evsi()
    #  if(is.null(input$Pop.OS)){return(NULL)}
    #  if(is.null(input$Time.OS)){return(NULL)}
    #  pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
    #  setup<-as.numeric(c(input$Setupmin,input$Setupmax))
    #  suppressWarnings(optimal<-optim.ss(evsi,setup,pp,input$Pop.OS,input$Time.OS,Dis=input$Dis,wtp=as.numeric(input$wtp.OS)))
    #  optimal$SS.max
    #  optimal$SS.I[2]
    #}
    #)
    output$ENBS<-shiny::renderText({
      evsi<-evsi()
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))

      
      paste("At the optimal sample size the Expected Net Benefit of Sampling is equal to ",
            round(suppressWarnings(optim.ss(evsi,setup,pp,as.numeric(input$Pop.OS),as.numeric(input$Time.OS),Dis=input$Dis,
                                            wtp=as.numeric(input$wtp.OS))$ENBS),-1),
      ". If this is greater than 0 then the study has economic benefit, if not then the ENBS demonstrates that the study not cost-effective.",sep="")
    })

    output$ENBS.plot<-shiny::renderPlot({
      evsi<-evsi()
      if(is.null(input$Pop.OS)){return(NULL)}
      if(is.null(input$Time.OS)){return(NULL)}
      pp<-as.numeric(c(input$PerPersmin,input$PerPersmax))
      setup<-as.numeric(c(input$Setupmin,input$Setupmax))

      enbs.plot(evsi,setup,pp,Pop=as.numeric(input$Pop.OS),Time=as.numeric(input$Time.OS),
                Dis=input$Dis,wtp=as.numeric(input$wtp.OS))
    })

    #output$Nmin<-shiny::renderText({
    #  evsi<-evsi()
    #  min(evsi$attrib$N)
    #})

    output$Nmax<-shiny::renderText({
      evsi<-evsi()
      paste("Note that the optimal sample size can only be found between",min(evsi$attrib$N),"and",
            max(evsi$attrib$N),
            "as these are the boundaries within which the EVSI has been calculated. If the optimal sample size is given as either of these values you will need to recalculate the EVSI for alternative values of N to find the true optimal sample size.",sep="")
    })


    #Reactive Values to give plots
    rv<-shiny::reactiveValues()
    shiny::observeEvent(input$Graphs,{rv$action<-"graphs"})
    shiny::observeEvent(input$Reload,{rv$action<-"reload"})

    output$main<-renderUI({
      panels<-list()
      data.panel<-{
        shiny::tabPanel("Data Upload",
                        shiny::tabsetPanel(
                          shiny::tabPanel("Model Upload",
                                          shiny::sidebarPanel(shiny::p("In this tab, you must upload the results of the baseline economic analysis.
                                                                       Firstly, the PSA simulations for the costs and effects must be uploaded. The simulations
                                                                       should be contained in the COLUMNS of a .csv file. Where the first two columns contain
                                                                       QALYs and then the costs for the first treatment and QALYS and costs for the next treatments
                                                                       are stored in the same manner in the next columns. An example of a spreadsheet in the correct
                                                                       form is given",shiny::a("here",href="http:://www.statistica.it/gianluca/BCEA/Vaccine_spreadsheet.csv")),
                                                              shiny::p("You must then upload the PSA simulations for the model parameters as a .csv file. The name
                                                                       of each parameter should be given as the first row of the spreadsheet."),
                                                              width=4),
                                          shiny::mainPanel(shiny::p(shiny::h4("Upload your files")),
                                                           shiny::fileInput("inputs","Upload a .csv file containing the PSA simulations for model parameters",
                                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                           shiny::fileInput("ec","Upload a .csv file containing the PSA simulations for costs and effects",
                                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                           width=8)),
                          shiny::tabPanel("EVSI Upload",
                                          shiny::sidebarPanel(
                                            shiny::p("In this tab you should upload a .csv files containing your calculated EVSI values.
                                                     The COLUMNS of the file should contain the EVSI for different willingness to pay values,
                                                     if the EVSI has only been calculated for one willingness to pay threshold then the
                                                     spreadsheet should only have one column. For ease, it is advisable to give the heading
                                                     of each column as the willingness to pay value for the EVSI."),
                                            shiny::p("The ROWS of the file should contain the EVSI for different sample sizes of the trial
                                                     under consideration. If the EVSI has only been calculated for one sample size then
                                                     the spreadsheet should only have one row. As before, it is advisable that the first
                                                     column contains the sample sizes for which the EVSI has been calculated."),
                                            width=4),
                                          shiny::mainPanel(
                                            shiny::p(shiny::h4("Upload your file")),
                                            shiny::fileInput("evsi.csv","Upload a .csv file containing the calculated EVSI values across willingness
                                                             to pay and sample size",
                                                             accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
                                            ,width=8)
                                            ),
                          shiny::tabPanel("EVSI Set Up",
                                          shiny::sidebarPanel(shiny::p("From the list below, please choose the parameters that your study is directly
                                                                       informing. For example, a clinical into the effectiveness of a novel treatment
                                                                       would directly inform the clinical outcome in the control and intervention arms"),
                                                              shiny::selectInput("parameter","Key parameters",choices=NULL,multiple=TRUE),
                                                              shiny::p("If the sample size for the study is NOT contained in the EVSI
                                                                       spreadsheet then please enter here, different sample sizes should be separated
                                                                       by a comma. Leave black is the sample size IS contained in the EVSI spreadsheet."),
                                                              shiny::textInput("N.input","Enter sample size",value=NULL),
                                                              shiny::p("Similarly, if the willingness to pay is NOT contained in the EVSI
                                                                       spreadsheet then please enter here, separated by a comma. Leave black is the
                                                                       sample size IS contained in the EVSI spreadsheet."),
                                                              shiny::textInput("wtp.input","Enter willingness to pay",value=NULL),
                                                              shiny::actionButton("RunModel","Run the analysis")
                                                              ,width=4),
                                          shiny::mainPanel(shiny::p("Once you have run the analysis, the standard cost-effectiveness analysis will be
                                                                    displayed here to check that the model has been uploaded correctly. For full
                                                                    economic analysis, use",shiny::a("BCEAweb.",href="https://egon.stats.ucl.ac.uk/projects/BCEAweb/")),
                                                           shiny::plotOutput("BCEA"),
                                                           shiny::actionButton("Graphs","Visualise the EVSI"),width=8)
                                          ,width=12)

                                                              ))
      }
      description.panel.data<-{
        #Description Panel
        shiny::tabPanel("Introduction to plotting",style='width:80%',
                        shiny::fluidRow(shiny::p("The Expected Value of Sample Information can be used to determine the economic benefit of a future
                                                 trial. However, in general, this economic benefit is conditional on several deterministic inputs that can
                                                 generally have a large impact on the value of the EVSI and any decision you make with the EVSI in mind.
                                                 Therefore, we have developed this tool to allow for the simple visualistion of the EVSI along with easy
                                                 manipulation of these deterministic inputs. This tool also allows for a simple user interface that can be
                                                 shown to collaborators and other stakeholders.",
                                                 shiny::h3("Willingness to Pay"),
                                                 shiny::p("The willingness to pay (WTP) is the amount of money the decision maker has avaliable to pay for an
                                                          additional 1 unit of benefit. This is typically given as a range and therefore the EVSI can be visualised
                                                          across different WTP values in the tab WTP. For the other analyses the WTP must be fixed. Although the plots
                                                          can easily be redrawn for different values of the WTP allowing the analysis to be easily completed for a
                                                          large number of different thresholds."),
                                                 shiny::h3("Sample Size"),
                                                 shiny::p("The sample size of the future trial is rarely determined before the EVSI analysis and one of the analysis
                                                          avaliable in this tool involves determining the optimal sample size for the future trial. The EVSI can be
                                                          visualised by sample size on the tab N. The value of a sample is bounded above by the value of resolving
                                                          all uncertainty in the parameters of interest. The key information to be gleaned from the by N plot is the
                                                          speed at which the EVSI reachs this upper bound."),
                                                 shiny::h3("Trial Cost-effectivenes"),
                                                 shiny::p("A key use of the EVSI is to determine whether a future trial will be cost-effective. This means that the
                                                          value of the future trial exceeds the cost of the trial. While in general this is a simple extension of the
                                                          ideas underpinning the EVSI, it typically depends on some additional inputs which are rarely known with
                                                          certainty.")),
                                        tags$ol(
                                          tags$li(shiny::p("Trial costs"),
                                                      shiny::p("The cost of the trial must be specified to allow comparison with the EVSI. However, these are rarely known
                                                               with certainty and so this application allows a range of values to be specified for the cost and then
                                                               considers the cost-effectiveness of the trial taking into account this uncertainty.")),
                                                      tags$li(shiny::p("Incidence Population"),
                                                      shiny::p("The EVSI is calculated as the EVSI per person who will benefit from the treatments under consideration.
                                                               Therefore, to compare with the trial costs, we must multiply by the number of people who will benefit
                                                               from the treatments. We call this the \"Incidence Population\". While there may be some literature that
                                                               can inform this parameter it cannot be known and therefore the cost-effectiveness of the trial is determined
                                                               for different levels of this population. This is an important consideration as the level of this population
                                                               can make a large difference to the cost-effectiveness of the trial.")),
                                                      tags$li(shiny::p("Time Horizon"),
                                                      shiny::p("Finally, the time horizon of a treatment is the length of time the treatments will be available in the
                                                               market before an alternative superseeds the treatments by being more cost-effective. Typically, this is
                                                               assumed to be around 10 years but we allow variation in this parameter to consider the cost-effectiveness
                                                               for alternative time horizons. This is powerful as it is rarely known the length of time a technology will
                                                               be available in the market."))),
                                                      actionButton("Reload","Reset Data")),
                        shiny::p("References: Raiffa, H. and Schlaifer, R.,",
                                 shiny::a("Applied Statistical Decision Theory, ",href="http://eu.wiley.com/WileyCDA/WileyTitle/productCd-047138349X.html"),
                                 shiny::em("Harvard University Press"),"1961"),
                        shiny::p("Wilson, E.,",
                                 shiny::a("A practical guide to value of information analysis, ",href="https://link.springer.com/article/10.1007%2Fs40273-014-0219-x"),
                                 shiny::em("PharmacoEconomics"),"2015"),
                        shiny::p("Steuten L., van de Wetering G., Groothuis-Oudshoorn K. and Retèl V.,",
                                 shiny::a("A Systematic and Critical Review of the Evolving Methods and Applications of Value of Information in Academia and Practice, ",
                                          href="https://link.springer.com/article/10.1007%2Fs40273-012-0008-3"),
                                 shiny::em("PharmacoEconomics"),"2013"),
                        shiny::p("The ",shiny::em("Visualisations for the EVSI "),"web app has been developed by ",shiny::a("Anna Heath", href="https://sites.google.com/site/annaheathstats/")," and ",
                                 shiny::a("Gianluca Baio",href="https://sites.google.com/a/statistica.it/gianluca/home"),
                                 "from the ",shiny::a("Statistics in Health Economic Evaluations Group",href="https://www.ucl.ac.uk/statistics/research/statistics-health-economics"),
                                 "at", shiny::a("University College London.",href="https://www.ucl.ac.uk/"),
                                 "Funding for the project was provided by the",shiny::a(" EPSRC",href="https://www.epsrc.ac.uk/"),
                                 "via a PhD Studentship and ",shiny::a("MAPI",href="http://mapigroup.com/"),
                                 "via an five year research grant."),
                        shiny::p("Copyright: ",shiny::a("Anna Heath", href="https://sites.google.com/site/annaheathstats/")," and ",
                                 shiny::a("Gianluca Baio",href="https://sites.google.com/a/statistica.it/gianluca/home"))
                        
                                                      )
      }
      description.panel<-{
        #Description Panel
        shiny::tabPanel("Introduction to plotting",style='width:80%',
                        shiny::fluidRow(shiny::p("The Expected Value of Sample Information can be used to determine the economic benefit of a future
                                                 trial. However, in general, this economic benefit is conditional on several deterministic inputs that can
                                                 generally have a large impact on the value of the EVSI and any decision you make with the EVSI in mind.
                                                 Therefore, we have developed this tool to allow for the simple visualistion of the EVSI along with easy
                                                 manipulation of these deterministic inputs. This tool also allows for a simple user interface that can be
                                                 shown to collaborators and other stakeholders.",
                                                 shiny::h3("Willingness to Pay"),
                                                 shiny::p("The willingness to pay (WTP) is the amount of money the decision maker has avaliable to pay for an
                                                          additional 1 unit of benefit. This is typically given as a range and therefore the EVSI can be visualised
                                                          across different WTP values in the tab WTP. For the other analyses the WTP must be fixed. Although the plots
                                                          can easily be redrawn for different values of the WTP allowing the analysis to be easily completed for a
                                                          large number of different thresholds."),
                                                 shiny::h3("Sample Size"),
                                                 shiny::p("The sample size of the future trial is rarely determined before the EVSI analysis and one of the analysis
                                                          avaliable in this tool involves determining the optimal sample size for the future trial. The EVSI can be
                                                          visualised by sample size on the tab N. The value of a sample is bounded above by the value of resolving
                                                          all uncertainty in the parameters of interest. The key information to be gleaned from the by N plot is the
                                                          speed at which the EVSI reachs this upper bound."),
                                                 shiny::h3("Trial Cost-effectivenes"),
                                                 shiny::p("A key use of the EVSI is to determine whether a future trial will be cost-effective. This means that the
                                                          value of the future trial exceeds the cost of the trial. While in general this is a simple extension of the
                                                          ideas underpinning the EVSI, it typically depends on some additional inputs which are rarely known with
                                                          certainty.")),
                                        tags$ol(
                                          tags$li(shiny::p("Trial costs"),
                                                  shiny::p("The cost of the trial must be specified to allow comparison with the EVSI. However, these are rarely known
                                                           with certainty and so this application allows a range of values to be specified for the cost and then
                                                           considers the cost-effectiveness of the trial taking into account this uncertainty.")),
                                          tags$li(shiny::p("Incidence Population"),
                                                  shiny::p("The EVSI is calculated as the EVSI per person who will benefit from the treatments under consideration.
                                                           Therefore, to compare with the trial costs, we must multiply by the number of people who will benefit
                                                           from the treatments. We call this the \"Incidence Population\". While there may be some literature that
                                                           can inform this parameter it cannot be known and therefore the cost-effectiveness of the trial is determined
                                                           for different levels of this population. This is an important consideration as the level of this population
                                                           can make a large difference to the cost-effectiveness of the trial.")),
                                          tags$li(shiny::p("Time Horizon"),
                                                  shiny::p("Finally, the time horizon of a treatment is the length of time the treatments will be available in the
                                                           market before an alternative superseeds the treatments by being more cost-effective. Typically, this is
                                                           assumed to be around 10 years but we allow variation in this parameter to consider the cost-effectiveness
                                                           for alternative time horizons. This is powerful as it is rarely known the length of time a technology will
                                                           be available in the market."))),
                                        shiny::p("References: Raiffa, H. and Schlaifer, R.,",
                                                 shiny::a("Applied Statistical Decision Theory, ",href="http://eu.wiley.com/WileyCDA/WileyTitle/productCd-047138349X.html"),
                                                 shiny::em("Harvard University Press"),"1961"),
                                        shiny::p("Wilson, E.,",
                                                 shiny::a("A practical guide to value of information analysis, ",href="https://link.springer.com/article/10.1007%2Fs40273-014-0219-x"),
                                                 shiny::em("PharmacoEconomics"),"2015"),
                                        shiny::p("Steuten L., van de Wetering G., Groothuis-Oudshoorn K. and Retèl V.,",
                                                 shiny::a("A Systematic and Critical Review of the Evolving Methods and Applications of Value of Information in Academia and Practice, ",
                                                          href="https://link.springer.com/article/10.1007%2Fs40273-012-0008-3"),
                                                 shiny::em("PharmacoEconomics"),"2013"),
                                        shiny::p("The ",shiny::em("Visualisations for the EVSI "),"web app has been developed by ",shiny::a("Anna Heath", href="https://sites.google.com/site/annaheathstats/")," and ",
                                                 shiny::a("Gianluca Baio",href="https://sites.google.com/a/statistica.it/gianluca/home"),
                                                 "from the ",shiny::a("Statistics in Health Economic Evaluations Group",href="https://www.ucl.ac.uk/statistics/research/statistics-health-economics"),
                                                 "at", shiny::a("University College London.",href="https://www.ucl.ac.uk/"),
                                                 "Funding for the project was provided by the",shiny::a(" EPSRC",href="https://www.epsrc.ac.uk/"),
                                                 "via a PhD Studentship and ",shiny::a("MAPI",href="http://mapigroup.com/"),
                                                 "via an five year research grant."),
                                        shiny::p("Copyright: ",shiny::a("Anna Heath", href="https://sites.google.com/site/annaheathstats/")," and ",
                                                 shiny::a("Gianluca Baio",href="https://sites.google.com/a/statistica.it/gianluca/home"))
                                                      )
        )
      }
      wtp.panel<-{
        shiny::tabPanel("EVSI by Willingness To Pay",style='width:80%',
                        shiny::sidebarPanel(shiny::p("The EVSI changes depending on the willingness-to-pay of the underlying decision maker and therefore
                                                     it is useful to visualise the EVSI for different WTP thresholds. This plot is also used visualise the
                                                     EVPI - which gives the maximum value for ANY future trial - and the EVPPI - which gives a maximum
                                                     value for a study targeting uncertainty in the parameters of interest of our study.
                                                     The EVSI also changes for different sample sizes as studies increase in value as the sample size increases.
                                                     The sample size can therefore be changed:"),
                                            shiny::selectInput(inputId="n",label="Choose a sample size",
                                                               choices=NULL,
                                                               selected=NULL),
                                            shiny::p("The EVSI should remain below the EVPPI for all sample sizes and willingness-to-pay thresholds.
                                                     All three measures should reach a sharp peak - this represents the \"break-even\" point between the
                                                     two treatment options where the decision uncertainty is at its maximum as the two treatment options are
                                                     equally likely to be cost-effective."),width=4),
                        shiny::mainPanel(shiny::plotOutput(outputId="ppEVSIbywtp"),
                                         shiny::p("References: Heath A. and Baio, G., ",
                                                  shiny::a("An Efficient Calculation Method for the Expected Value of Sample Information: Can we do it? Yes, we can ",
                                                           href="https://arxiv.org/abs/1709.02319"),
                                                  shiny::em(", arXiv preprint"),", 2018"),
                                         shiny::p("McCabe, C., Claxton K. and Culyer, A., ",
                                                  shiny::a("The NICE cost-effectiveness threshold",href="https://link.springer.com/article/10.2165/00019053-200826090-00004"),
                                                  shiny::em(", PharmacoEconomics"),", 2008"),
                                         shiny::p("Baio, G., ",
                                                  shiny::a("Bayesian Methods in Health Economics",href="https://sites.google.com/a/statistica.it/gianluca/bookhe"),
                                                  shiny::em(", Springer"),", 2012"),
                                         shiny::p("Baio, G., Berardi, A. and Heath, A., ",
                                                  shiny::a("Bayesian Cost-Effectiveness Analysis with the R package BCEA",href="http://www.springer.com/gb/book/9783319557168"),
                                                  shiny::em(", Springer"),", 2017"),width=8)
                                            )

      }
      N.panel<-{
        shiny::tabPanel("EVSI by Sample Size",style='width:80%',
                        shiny::sidebarPanel(shiny::p("The EVSI increases as the sample size of the underlying trial increases. This graphic shows the
                                                     EVSI across different sample sizes. This relationship with N changes depending on the value of the
                                                     willingness-to-pay so the plot can be considered for changing values of the WTP.
                                                     The EVSI is calculated using Bayesian regression. This means that posterior credible intervals for the EVSI
                                                     can be calculated and these are plotted on the graphic. This demonstrates the uncertainty in the EVSI
                                                     estimate. If the uncertainty is too large for decision making then Q should be increased in the
                                                     EVSI calculation."),
                                            shiny::selectInput(inputId="wtp",label="Choose a Willingness-to-Pay Threshold",
                                                               choices=NULL,selected=NULL),
                                            shiny::p(paste("In general, the EVSI will be more accurately estimated for higher values of the EVSI. This
                                                           would will be for willingness to pay values close to the \"break-even\" point.")),width=4),
                        shiny::mainPanel(shiny::plotOutput(outputId="ppEVSIbyn"),
                                         shiny::p("References: Heath, A., Manolopoulou I. and Baio, G., ",
                                                  shiny::a("Efficient Monte Carlo Estimation of the Expected Value of Sample Information using Moment Matching",href="https://arxiv.org/abs/1611.01373"),
                                                  shiny::em(", Medical Decision Making"),", forthcoming"),
                                         shiny::p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                  shiny::a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                  shiny::em(", Working Paper"),", 2017"),
                                         width=8)
                        )
      }
      CE.panel.multi.N<-{
        shiny::tabPanel("Cost-effectiveness of a Trial",style='width:80%',
            shiny::column(12,
                  shiny::tabsetPanel(id="CE",
                        #Set up Panel
                        shiny::tabPanel("Setup",
                                  shiny::fluidRow(
                                     shiny::sidebarPanel({
                                       shiny::p("A trial is const-effective if the value of information gained from the trial exceeds the cost of undertaking the trial.
                                                Therefore, to determine whether the proposed trial is cost-effective we must first determine the cost of the trial under consideration.
                                                This tab allows you to specify key characteristics of the trial so the cost-effectiveness of the trial can be ascertained.
                                                The following tabs display the cost-effectiveness of the trial graphically and, if the EVSI has been calculated for different sample sizes,
                                                allows you to determine the optimal sample size of the trial. For more information about decision making using the EVSI in clinical trials,
                                                see the references below.")}),
                                          shiny::mainPanel({
                                             shiny::tabsetPanel(id="set",
                                                 shiny::tabPanel("Trial Costs",
                                                                  shiny::fluidRow(shiny::p("In general, trial costs are split into two categories; the setup and the per person costs.
                                                                                                 The are setup costs are the overhead costs of the trial and will be incurred irrespective of the size
                                                                                                 of the trial. Typical setup costs may be training for staff or the purchase of specialised equipment.
                                                                                                 Additional costs will then be incurred for each participant enrolled in the trial such as the cost of
                                                                                                 administering the treatment or following up the patient."),
                                                                          shiny::p("In most settings, these two costs will not be known with certainty. Therefore, you should give a range of possible
                                                                               values for these costs. If the costs are known with certainty then simply input the same values for
                                                                                   the maximum and minimum possible values of the costs. ")),
                                                                          shiny::fluidRow(
                                                                            shiny::column(4,shiny::numericInput(inputId="Setupmin",label="Minimum Setup Costs for the Trial",
                                                                                                                value=150000,step=10,min=0),
                                                                                          shiny::numericInput(inputId="PerPersmin",label="Minimum Cost Per Person",min=0,
                                                                                                              value=500,step=10)),
                                                                            shiny::column(4,shiny::numericInput(inputId="Setupmax",label="Maximum Setup Costs for the Trial",
                                                                                                                value=1500000,step=10,min=0),
                                                                                          shiny::numericInput(inputId="PerPersmax",label="Maximum Cost Per Person",min=0,
                                                                                                              value=2000,step=10)
                                                                            )))
                                                                        ,
                                                 shiny::tabPanel("Long Term Dynamics",
                                                     shiny::fluidRow(shiny::p("The cost-effectiveness of a trial depends on two additional inputs. These are the incidence population, i.e.
                                                                              the yearly incidence of the disease under consideration. This gives the number of patients that will benefit from
                                                                                                               the treatment in each year that the treatment is used. This may not be known with certainty and the plot considers
                                                                                                               the cost-effectiveness of the trial for different possible numbers of patients. However, it is necessary to give
                                                                                                               possible values of the incidence population."),
                                                                     shiny::p("The time horizon gives the number of years that before a more effective/efficient treatment will enter the market.
                                                                                                               This can be thought of as the number of years before a more effective treatment will be developed.
                                                                                                               This will depend on the disease areas as fast moving diseases such as cancer will have a shorter time
                                                                                                               horizon. The maximum and minimum possible values for the time horizon should be specified here.")),
                                                        shiny::fluidRow(shiny::column(6,shiny::numericInput(inputId="Popmin",label="Minimum Incidence Population",value=0,step=100,min=0),
                                                                        shiny::numericInput(inputId="Timemin",label="Minimum Time Horizon",value=0,step=1,min=0)),
                                                                        shiny::column(6,
                                                                               shiny::numericInput(inputId="Popmax",label="Maximum Incidence Population",value=1e+05,step=100,min=0),
                                                                                shiny::numericInput(inputId="Timemax",label="Maximum Time Horizon",min=0,value=10,step=1)
                                                                                                      )),
                                                                                      shiny::fluidRow(shiny::p("Additional References: Thokala P., Goodacre S., Ward M., Penn-Ashman J. and Perkins G., ",
                                                                                                              shiny::a("Cost-effectiveness of Out-of-Hospital Continuous Positive Airway Pressure for Acute Respiratory Failure.",
                                                                                                                       href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4414542/"),
                                                                                                              shiny::em(", Annals of Emergency Medicine"),", 2015"),
                                                                                                      shiny::p("McKenna C. and Claxton K., ",
                                                                                                               shiny::a(" Addressing Adoption and Research Design Decisions Simultaneously: The Role of
                                                                                                                        Value of Sample Information Analysis",href="http://journals.sagepub.com/doi/pdf/10.1177/0272989X11399921"),
                                                                                                               shiny::em(", Medical Decision Making"),", 2011"))
                                                                                      ),
                                                 shiny::tabPanel("Discount Rate",
                                                          shiny::fluidRow(shiny::p("The final input to determine the cost-effectiveness of the trial is the discount rate for future treatments.
                                                                                                 In general, health benefits now are more valuable than health benefits in the future. NICE recommend 3.5% as the
                                                                                               discount rate for the treatments but this can be changed here."),
                                                           shiny::numericInput(inputId="Dis",label="Discount Rate",
                                                                                                          value=0.035,step=0.001),
                                                           shiny::p("Additional References: NICE,",
                                                                    shiny::a(" Methods for the development of NICE public health guidance (third edition)",
                                                                             href="https://www.nice.org.uk/process/pmg4/chapter/incorporating-health-economics"),
                                                                    shiny::em("NICE Guidelines"),", 2012"))
                                                                                      )
                                                 )})
                                                            ),
                                   shiny::fluidRow({shiny::column(8,offset=4,shiny::p("References: Willan A. and Pinto E., ",
                                                                                     shiny::a("The value of information and optimal clinical trial design.",
                                                                                        href="http://onlinelibrary.wiley.com/doi/10.1002/sim.2069/pdf"),
                                                                                     shiny::em(", Statistics in Medicine"),", 2005"),
                                                                                     shiny::p("Eckermann S. and Willan A., ",
                                                                                             shiny::a(" The option value of delay in health technology assessment",href="http://journals.sagepub.com/doi/pdf/10.1177/0272989X07312477"),
                                                                                             shiny::em(", Medical Decision Making"),", 2008") 
                                                                                                      )})
                        ),
                        #Cost-effective Trial
                        shiny::tabPanel("Probability of CE Trial",
                                                                                        shiny::fluidRow(shiny::sidebarPanel(#Population EVSI - prob of CE plot
                                                                                          shiny::selectInput(inputId="n.CE",label="Choose a sample size",
                                                                                                             choices=NULL,selected=NULL),
                                                                                          shiny::selectInput(inputId="wtp.CE",label="Choose a Willingness-to-Pay Threshold",
                                                                                                             choices=NULL,selected=NULL),
                                                                                          shiny::uiOutput("PopDynam"),
                                                                                          shiny::uiOutput("TimeDynam"),width=4),
                                                                                          shiny::mainPanel(shiny::plotOutput(outputId="ProbCE"),
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
                                                                                                                      shiny::p("References: Heath, A., Baio, G., and Hunter, R. ",
                                                                                                                               shiny::a("Development of a New Software Tool to Compute the Expected Value of Sample Information - An application to the HomeHealth intervention",href="https://sites.google.com/site/annaheathstats"),
                                                                                                                               shiny::em(", HESG Winter Meeting"),", 2018"),
                                                                                                                      shiny::p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                                               shiny::a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                                                                                               shiny::em(", Working Paper"),", 2017")))
                                                                                        ),
                        #Optimal Sample Size
                        shiny::tabPanel("Optimal Sample Size",value = "OSS",
                                                                                        shiny::sidebarPanel(shiny::selectInput(inputId="wtp.OS",label="Choose a Willingness-to-Pay Threshold",
                                                                                                                               choices=NULL,selected=NULL),
                                                                                                            shiny::uiOutput("Pop.OSDynam"),
                                                                                                            shiny::uiOutput("Time.OSDynam"),
                                                                                                            shiny::p("It is possible to find the sample size for your trial that will give the maximum value for money."),
                                                                                                            shiny::p(shiny::textOutput(outputId = "Nmax")),
                                                                                                            shiny::p("The plot for the ENBS allows you to assess how the ENBS is behaving for different sample sizes. Please
                                                                                                                     be aware that it is often possible for a large number of sample sizes to have similar economic value and therefore
                                                                                                                     the optimal sample size should be interpreted with care.")
                                                                                                            ,width=4),

                                                                                        shiny::mainPanel(shiny::fluidRow(shiny::plotOutput("ENBS.plot"),
                                                                                                                         shiny::column(5,shiny::p(shiny::textOutput(outputId="SS"))),
                                                                                                                         shiny::column(7,shiny::p(shiny::textOutput(outputId="ENBS"))),
                                                                                                                         width=8))
                                                                                                                         )
          ))
        )
        
      }
      CE.panel.single.N<-{
        shiny::tabPanel("Cost-effectiveness of a Trial",style='width:80%',
                        shiny::column(12,
                              shiny::tabsetPanel(id="CE",
                                                 shiny::tabPanel("Setup",
                                                                 shiny::fluidRow(
                                                                   shiny::sidebarPanel({
                                                                     shiny::p("A trial is const-effective if the value of information gained from the trial exceeds the cost of undertaking the trial.
                                                                              Therefore, to determine whether the proposed trial is cost-effective we must first determine the cost of the trial under consideration.
                                                                              This tab allows you to specify key characteristics of the trial so the cost-effectiveness of the trial can be ascertained.
                                                                              The following tabs display the cost-effectiveness of the trial graphically and, if the EVSI has been calculated for different sample sizes,
                                                                              allows you to determine the optimal sample size of the trial. For more information about decision making using the EVSI in clinical trials,
                                                                              see the references below.")}),
                                                                   shiny::mainPanel({
                                                                     shiny::tabsetPanel(id="set",
                                                                                        shiny::tabPanel("Trial Costs",
                                                                                                        shiny::fluidRow(shiny::p("In general, trial costs are split into two categories; the setup and the per person costs.
                                                                                                                                 The are setup costs are the overhead costs of the trial and will be incurred irrespective of the size
                                                                                                                                 of the trial. Typical setup costs may be training for staff or the purchase of specialised equipment.
                                                                                                                                 Additional costs will then be incurred for each participant enrolled in the trial such as the cost of
                                                                                                                                 administering the treatment or following up the patient."),
                                                                                                                        shiny::p("In most settings, these two costs will not be known with certainty. Therefore, you should give a range of possible
                                                                                                                                 values for these costs. If the costs are known with certainty then simply input the same values for
                                                                                                                                 the maximum and minimum possible values of the costs. ")),
                                                                                                        shiny::fluidRow(
                                                                                                          shiny::column(4,shiny::numericInput(inputId="Setupmin",label="Minimum Setup Costs for the Trial",
                                                                                                                                              value=150000,step=10,min=0),
                                                                                                                        shiny::numericInput(inputId="PerPersmin",label="Minimum Cost Per Person",min=0,
                                                                                                                                            value=500,step=10)),
                                                                                                          shiny::column(4,shiny::numericInput(inputId="Setupmax",label="Maximum Setup Costs for the Trial",
                                                                                                                                              value=1500000,step=10,min=0),
                                                                                                                        shiny::numericInput(inputId="PerPersmax",label="Maximum Cost Per Person",min=0,
                                                                                                                                            value=2000,step=10)
                                                                                                          )))
                                                                                        ,
                                                                                        shiny::tabPanel("Long Term Dynamics",
                                                                                                        shiny::fluidRow(shiny::p("The cost-effectiveness of a trial depends on two additional inputs. These are the incidence population, i.e.
                                                                                                                                 the yearly incidence of the disease under consideration. This gives the number of patients that will benefit from
                                                                                                                                 the treatment in each year that the treatment is used. This may not be known with certainty and the plot considers
                                                                                                                                 the cost-effectiveness of the trial for different possible numbers of patients. However, it is necessary to give
                                                                                                                                 possible values of the incidence population."),
                                                                                                                        shiny::p("The time horizon gives the number of years that before a more effective/efficient treatment will enter the market.
                                                                                                                                 This can be thought of as the number of years before a more effective treatment will be developed.
                                                                                                                                 This will depend on the disease areas as fast moving diseases such as cancer will have a shorter time
                                                                                                                                 horizon. The maximum and minimum possible values for the time horizon should be specified here.")),
                                                                                                        shiny::fluidRow(shiny::column(6,shiny::numericInput(inputId="Popmin",label="Minimum Incidence Population",value=0,step=100,min=0),
                                                                                                                                      shiny::numericInput(inputId="Timemin",label="Minimum Time Horizon",value=0,step=1,min=0)),
                                                                                                                        shiny::column(6,
                                                                                                                                      shiny::numericInput(inputId="Popmax",label="Maximum Incidence Population",value=1e+05,step=100,min=0),
                                                                                                                                      shiny::numericInput(inputId="Timemax",label="Maximum Time Horizon",min=0,value=10,step=1)
                                                                                                                        )),
                                                                                                        shiny::fluidRow(shiny::p("Additional References: Thokala P., Goodacre S., Ward M., Penn-Ashman J. and Perkins G., ",
                                                                                                                                 shiny::a("Cost-effectiveness of Out-of-Hospital Continuous Positive Airway Pressure for Acute Respiratory Failure.",
                                                                                                                                          href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4414542/"),
                                                                                                                                 shiny::em(", Annals of Emergency Medicine"),", 2015"),
                                                                                                                        shiny::p("McKenna C. and Claxton K., ",
                                                                                                                                 shiny::a(" Addressing Adoption and Research Design Decisions Simultaneously: The Role of
                                                                                                                                          Value of Sample Information Analysis",href="http://journals.sagepub.com/doi/pdf/10.1177/0272989X11399921"),
                                                                                                                                 shiny::em(", Medical Decision Making"),", 2011"))
                                                                                                        ),
                                                                                        shiny::tabPanel("Discount Rate",
                                                                                                        shiny::fluidRow(shiny::p("The final input to determine the cost-effectiveness of the trial is the discount rate for future treatments.
                                                                                                                                 In general, health benefits now are more valuable than health benefits in the future. NICE recommend 3.5% as the
                                                                                                                                 discount rate for the treatments but this can be changed here."),
                                                                                                                        shiny::numericInput(inputId="Dis",label="Discount Rate",
                                                                                                                                            value=0.035,step=0.001),
                                                                                                                        shiny::p("Additional References: NICE,",
                                                                                                                                 shiny::a(" Methods for the development of NICE public health guidance (third edition)",
                                                                                                                                          href="https://www.nice.org.uk/process/pmg4/chapter/incorporating-health-economics"),
                                                                                                                                 shiny::em("NICE Guidelines"),", 2012"))
                                                                                                        )
                                                                                                        )})
                                                                                                        ),
                                                                 shiny::fluidRow({shiny::column(8,offset=4,shiny::p("References: Willan A. and Pinto E., ",
                                                                                                                    shiny::a("The value of information and optimal clinical trial design.",
                                                                                                                             href="http://onlinelibrary.wiley.com/doi/10.1002/sim.2069/pdf"),
                                                                                                                    shiny::em(", Statistics in Medicine"),", 2005"),
                                                                                                shiny::p("Eckermann S. and Willan A., ",
                                                                                                         shiny::a(" The option value of delay in health technology assessment",href="http://journals.sagepub.com/doi/pdf/10.1177/0272989X07312477"),
                                                                                                         shiny::em(", Medical Decision Making"),", 2008") 
                                                                 )})
                                                                   ),
                                                                                         shiny::tabPanel("Probability of CE Trial",
                                                                                                         shiny::fluidRow(shiny::sidebarPanel(#Population EVSI - prob of CE plot
                                                                                                           shiny::selectInput(inputId="n.CE",label="Choose a sample size",
                                                                                                                              choices=NULL,selected=NULL),
                                                                                                           shiny::selectInput(inputId="wtp.CE",label="Choose a Willingness-to-Pay Threshold",
                                                                                                                              choices=NULL,selected=NULL),
                                                                                                           shiny::uiOutput("PopDynam"),
                                                                                                           shiny::uiOutput("TimeDynam"),width=4),
                                                                                                           shiny::mainPanel(shiny::plotOutput(outputId="ProbCE"),
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
                                                                                                                                       shiny::p("References: Heath, A., Baio, G., and Hunter, R.",
                                                                                                                                                shiny::a("Development of a New Software Tool to Compute the Expected Value of Sample Information - An application to the HomeHealth intervention",href="https://sites.google.com/site/annaheathstats"),
                                                                                                                                                shiny::em(", HESG Winter Meeting"),", 2018"),
                                                                                                                                       shiny::p("Heath, A., Manolopoulou I. and Baio, G., ",
                                                                                                                                                shiny::a(" Bayesian Curve Fitting to Estimate the Expected Value of Sample Information using Moment Matching Across Different Sample Sizes",href="https://sites.google.com/site/annaheathstats/selected-publications/curve-fitting-paper"),
                                                                                                                                                shiny::em(", Working Paper"),", 2017")))
                                                                                                         )

                                                                                                           ))
                                                                                                         )


      }

      if(complex=="y"){panels[[1]]<-data.panel

      if(length(rv$action=="graphs")>0){
        if(rv$action=="graphs"){
          panels[[1]]<-description.panel.data
          panels[[2]]<-wtp.panel

          if(length(evsi()$attrib$N)>1){
            panels[[3]]<-N.panel
            panels[[4]]<-CE.panel.multi.N
          }
          if(length(evsi()$attrib$N)==1){
            panels[[3]]<-CE.panel.single.N
          }
        }
        if(rv$action=="reload"){
          panels<-list(data.panel)
        }
      }
      }
      if(simple=="y"){
        panels[[1]]<-description.panel
        panels[[2]]<-wtp.panel

        if(length(evsi()$attrib$N)>1){
          panels[[3]]<-N.panel
          panels[[4]]<-CE.panel.multi.N
        }
        if(length(evsi()$attrib$N)==1){
          panels[[3]]<-CE.panel.single.N
        }
      }
      do.call(shiny::tabsetPanel,panels)
      })


  }

  suppressWarnings(shiny::shinyApp(ui=ui,server=server))

    }

