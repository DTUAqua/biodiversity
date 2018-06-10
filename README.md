# biodiversity
Git page for `biodiversity` R-package containing models for biodiversity 

**Running the GAM analysis
 
The source code for the GAM analysis is available in the file in the github repository called:
 
 extra.R
 
To run the analysis in this file require that R (http://r-project.org) and the supporting R-packages are installed. The needed packages are listed in the file extra.R in the lines 4-16. Each of these can be installed via the install.packages function in R. E.g. the needed package "gtable" can be installed via the command: 
 
 install.packages("gtable")
 
After all required packages have been installed it will be possible to replicate the analysis by downloading the file extra.R and then sourcing the file in R via the command: 
 
 source("extra.R")    
 
 
**Installing the biodiversity package
 
The biodiversity models implemented for this manuscript have been collected in an R-package named "biodiversity". This package is available via the github site. The package depends on an open source and freely available package for automatic differentiation (http://tmb-project.org). Instructions for installing this package are collected here: 
 
 https://github.com/kaskr/adcomp/wiki/Download  
 
To install the package from the password protected github site require the devtools package. This can be installed in R via the command: 
 
 install.packages("devtools")
 
After this the biodiversity package can be installed via the line: 
 
 devtools::install_github("DTUAqua/biodiversity/biodiversity", auth_token="41e59f73cfe7b855dd61b39db125b79cb235d9d7")
 
The package contain all needed data and functions to run the different models - these are documented in the package. E.g. to run the metabolic model type in R: 
 
 library("biodiversity")
 data(species)
 fit <- biodiv(species, conf=4)
 plot(fit)
 coef(fit)
 
Here the last line will result in a table of the model parameter estimates. 
 
The authors of the package are available to answer any questions the reviewers may have. 

