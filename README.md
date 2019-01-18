[![Travis-CI Build Status](https://travis-ci.org/youjin1207/netchain.svg?branch=master)](https://travis-ci.org/youjin1207/netchain)
![![Downloads badge](http://cranlogs.r-pkg.org/badges/netchain)](http://cranlogs.r-pkg.org/badges/netchain?color=red)
 [![](http://cranlogs.r-pkg.org/badges/grand-total/netchain?color=yellow)](https://CRAN.R-project.org/package=netchain)
[![arXiv shield](https://img.shields.io/badge/arXiv-1812.04990-blue.svg?style=flat)](https://arxiv.org/abs/1812.04990)

# Overview

netchain is a R package for causal inference on collective outcomes under social network. [Our paper](https://arxiv.org/abs/1812.04990) proposed and justified a parsimonious parametrization for social network data generated from causal directed acyclic graph (DAG), approximating a particular family of graphical models known as chain graphs under some conditions. 

We provide a function `simGibbs()` to generate binary outcomes, treatments, and confounders from chain graph model. A function `chain.causal.multi()` is to infer parameters in the conditional log-linear models that feature hybrid graphical models of undirected graphs and directed acyclic graphs (DAG). This function generates counterfactual outcomes using Gibbs sampling given `treatment` assignment and the estimated parameters to derive the probability associated with collective outcomes. We also provide a function of `causal.influence()` to identify the most (causally) influential subjects in social network based on the their causal effect on the collective outcomes. 

## Package information

- Version: 0.1.0
- Author : Elizabeth Ogburn (<eogburn@jhu.edu>), Ilya Shpitser (<ilyas@jhu.edu>), and Youjin Lee (<ylee160@jhu.edu>)
- Maintainer : Youjin Lee (<ylee160@jhu.edu>)
- Imports: Rcpp (>= 0.12.17), Matrix, gtools, stringr, stats, igraph, simcausal
- Linking to : Rcpp

## Installation

You can download the package by:
```
install.packages("netchain")

# or you can directly download the development version from author's Github
install.packages("devtools")
library(devtools)
install_github("youjin1207/netchain")
```

## Usage


[Here](https://github.com/youjin1207/netchain/blob/master/vignettes/chainapprox.Rmd) is a R vignettes for guidance. Or you can access to vignettes via:

```
install_github("youjin1207/netchain", build_vignettes = TRUE)
vignette("chainapprox", package = "netchain")
```

### Example

```
library(netchain)
# set direct effect and two-way interaction effect on undirected graphs (weight.matrix)
weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
simobs = simGibbs(n.unit = 3, n.gibbs = 10, n.sample = 10, 
                   weight.matrix,
                   treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
inputY = simobs$inputY
inputA = simobs$inputA
inputC = simobs$inputC

# define relational matrix (R.matrix)
R.matrix = ifelse(weight.matrix==0, 0, 1)      
diag(R.matrix) = 0

# infer conditional log-linear model following chain graph models.
result = chain.causal.multi(targetoutcome = "mean", treatment = c(1,0,0), inputY, inputA, listC = inputC, R.matrix = R.matrix, E.matrix = diag(3), edgeinfo = list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), rbind(c("Y", 3), c("C", 3))), n.obs = 1000, n.burn = 100)
print(result)

# measure influence for each node by evaluating average of collective outcomes under each treatment.
influence = causal.influence(targetoutcome = "mean", Avalues = c(1,0), 
                            inputY, inputA, listC = inputC, R.matrix, E.matrix = diag(3), 
                            edgeinfo = list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), rbind(c("Y", 3), c("C", 3))), n.obs = 100, n.burn = 10)
print(influence)
```

## Reference

Ogburn, E. L., Shpitser, I., \& Lee, Y. (2018). Causal inference, social networks, and chain graphs. _arXiv preprint arXiv:1812.04990_.
