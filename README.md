# EVSI
Based on a Bayesian health economic model which includes the current information about the model parameters and distributions for the future data collection, functions are available to calculate the Expected Value of Sample Information (EVSI) using the Heath et al. Moment Matching method (2017).
This package then contains a number of plots that present the EVSI and the Expected Net Benefit of Sampling (ENBS) to present and analysis the EVSI calculation. Finally, a `shiny` web application is available to dynamically explore the EVSI and present results to key stakeholders.

# Installation
The `EVSI` package is currently only available in a developement version through GitHub and must be installed using the package `devtools`
```R
install.packages("devtools")
```
The EVSI calculation method is based on Bayesian updating so an MCMC sampler is needed. Currently the EVSI package can be used with either [jags](http://mcmc-jags.sourceforge.net/) or [OpenBUGS](http://www.openbugs.net/w/FrontPage). These need to be installed separately from their respective repositories and instuctions for installations under different OS can be found online.
Depending on which MCMC software is used the EVSI package either requires the package `rjags` or the pacakge `R2OpenBUGS`, note that both these packages require that the respective MCMC sampler is installed separately
```R
install.packages("rjags")
install.packages("R2OpenBUGS")
```
In addition to the MCMC sampler the following packages must be installed before proceeding to EVSI calculation and presentation 
```R
install.packages("BCEA","shiny","shinythemes")
```

After installing the required dependancies, the `EVSI` package can be installed using `devtools` with the following code
```R
devtools::install_github("annaheath/EVSI")
```
Note that the `EVSI` package is currently under active developement and therefore it is advisable to reinstall the package directly from GitHub before each use to ensure that you are using the most uptodate version.
