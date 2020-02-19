# biodiversity
Git page for `biodiversity` R-package containing models for biodiversity 

### Installing the biodiversity package
 
The biodiversity models implemented for this manuscript have been collected in an R-package named "biodiversity". This package is available via the github site. The package depends on an open source and freely available package for automatic differentiation (http://tmb-project.org). Instructions for installing this package are collected here: 
 
 https://github.com/kaskr/adcomp/wiki/Download  
 
To install the package from the password protected github site require the devtools package. This can be installed in R via the command: 

```R
 install.packages("devtools")
```

After this the biodiversity package can be installed via the line: 

```R
 devtools::install_github("DTUAqua/biodiversity/biodiversity")
```

The package contain all needed data and functions to run the different models - these are documented in the package. E.g. to run the metabolic model type in R: 

```R
 library("biodiversity")
 data(species)
 fit <- biodiv(species, conf=4)
 plot(fit)
 coef(fit)
```

Here the last line will result in a table of the model parameter estimates. 

The results in table 1 of Gislason et al. (2020) can be easily reproduced via the script: 
```R
test/loglik/script.R
```

# Reference

Gislason et al. 2020. Species richness in North Atlantic fish: Process concealed by pattern. Global Ecology and Biogeography. DOI: 10.1111/geb.13068
