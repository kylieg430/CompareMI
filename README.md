
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CompareMI

<!-- badges: start -->

[![R-CMD-check](https://github.com/kylieg430/CompareMI/workflows/R-CMD-check/badge.svg)](https://github.com/kylieg430/CompareMI/actions)
<!-- badges: end -->

The goal of CompareMI is to compare multiple imputation methods using
RMSE and clinically relevant measures.

## Installation

You can install the released version of CompareMI from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CompareMI")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kylieg430/CompareMI")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CompareMI)
Data <- generateData(100)
path <-  "C:/Users/kgetz1/AppData/Local/Programs/Python/Python38/python.exe"
workingdir <-  "C:/Users/kgetz1/Documents/Year2/ComputingProject/GitHub/CompareMI"
runSims("MCAR",size=100,iterations=2,sdat.c=Data,parallelOption=FALSE,pythonPath=path,wd=workingdir)
runSims("MCAR",size=100,iterations=2,sdat.c=Data,parallelOption=TRUE,coreNum=2,pythonPath=path,wd=workingdir)
RestructureResults(mech="MCAR",size=100,wd=workingdir)
plotResults("MCAR",-4,4,size=100,wd=workingdir)
calcRMSE("MCAR",size=100,sdat.c=Data,pythonPath=path,wd=workingdir)
```
