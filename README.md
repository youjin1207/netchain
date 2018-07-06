[![Travis-CI Build Status](https://travis-ci.org/youjin1207/netchain.svg?branch=master)](https://travis-ci.org/youjin1207/netchain)

# netchain
R package for causal inference on collective outcomes under social network

- Version: 0.1.0
- Author : Youjin Lee (<ylee160@jhu.edu>)
- Maintainer : Youjin Lee (<ylee160@jhu.edu>)
- Imports : Rcpp (>= 0.12.17), Matrix, gtools, stringr

You can download the package by:
```
install.packages("devtools")
library(devtools)
install_github("youjin1207/netchain")
```
[Here](https://github.com/youjin1207/netchain/blob/master/vignettes/chainapprox.Rmd) is a R vignettes for guidance. Or you can access to vignettes via:

```
install_github("youjin1207/netchain", build_vignettes = TRUE)
vignette("chainapprox", package = "netchain")
```
