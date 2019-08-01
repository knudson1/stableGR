# stableGR


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)

# Overview
Welcome to the stableGR repository, containing the production-ready version of the R package `stableGR`.

# Repo Contents

- [R](./stableGR/R): `R` package code.
- [man](./stableGR/man): package documentation and usage of the `stableGR` package.
- [data](./stableGR/data): data included in the`stableGR` package.


# System Requirements

## Hardware Requirements

The `stableGR` package requires no special hardware.  Because the `stableGR` package diagnoses convergence for Markov chains, computers with greater RAM will complete the analyses more quickly.




## OS Requirements

The `stableGR` package was developed and tested on Ubuntu using version 3.4.3 of R. Additional testing was performed on a Windows PC and on a Mac. The package should be compatible with Windows, Mac, and Linux operating systems.









# Installation Guide

Before installing the `stableGR` package, users should have downloaded and installed `R` version 3.1.1 or higher from from [CRAN](https://cran.r-project.org/).  Installing R  should take 2 to 3 minutes and installing the `stableGR` package should take 1 to 2 minutes. 

This package will be available on [CRAN](https://cran.r-project.org/) once updates from the R package `mcmcse` are pushed to CRAN. In the meantime, the `devtools` package enables users to  install packages directly from GitHub. 

First, you will need to get three required packages:
```
install.packages("Rcpp") 
install.packages("RcppArmadillo")
install.packages("devtools")
```
Then, install mcmcse from github  (rather than CRAN) so that you have the newer methods, which are not yet available on CRAN.
```
library(devtools)
install_github("dvats/mcmcse")
```

Finally, you can install the `stableGR` package:
```
install_github("knudson1/stableGR/stableGR")
library(stableGR)
```

The above steps need to be performed once. Once `stableGR` is installed on your computer, you can simply type
```
library(stableGR)
```
to use the package (just like if you had installed a package from CRAN).

Detailed directions can be found online in many places, including [this source](http://kbroman.org/pkg_primer/pages/github.html).







# License

The license for `stableGR` is [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html).




