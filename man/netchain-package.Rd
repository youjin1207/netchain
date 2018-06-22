\name{netchain-package}
\alias{netchain-package}
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
\references{
  
}
\keyword{ package }
\seealso{
  \url{https://github.com/youjin1207/netchain}
}
\examples{
set.seed(123)
weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
simobs = simGibbs(n.unit = 3, n.obs = 200, 
                   weight.matrix,
                  treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
inputY = simobs$inputY                   
inputA = simobs$inputA   
inputC = simobs$inputC 
R.matrix = ifelse(weight.matrix==0, 0, 1)      
result = chain.causal.multi(targetoutcome = "mean", treatment = c(1,0,0), 
                            inputY, inputA, listC = inputC, R.matrix, E.matrix = diag(3), 
                            edgeinfo = list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), rbind(c("Y", 3), c("C", 3))), 
                            n.obs = 1000, n.burn = 100)

}
