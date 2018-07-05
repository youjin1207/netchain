\name{netchain-package}
\alias{netchain}
\docType{package}
\title{
  netchain: causal inference on collective outcomes 
}
\description{
 This package is for estimation of probability associated with collective counterfactual outcomes using approximation via causal graphical model. We apply a parsimonious parameterization for social network data with some specific kinds of interference and contagion, which corresponds to particular family of graphical models known as chain graphs. 
}
\details{
  We provide functions to estimate the parameters in conditional log-linear model when the observations (outcomes, treatments, and confounders) and the structure of a causal graph is given. Based on the estimated parameters, we generate counterfactual outcomes using Gibbs sampling to infer the causal effect (or causal probability) of a certain treatment assignment on the collective outcomes. Moreover, we use this method to identify causally influential units on social network. 
}
\author{
Youjin Lee

Maintainer: Youjin Lee <ylee160@jhu.edu>
}

\keyword{ package }
\seealso{
  \url{https://github.com/youjin1207/netchain}
}

