#' Generate Gibbs samplers for counterfactual collective outcomes.
#'
#' This function generates the outcomes using Gibbs sampling under the given treatment assignment
#' and edge information. 
#'
#' @param pars a set of parameters
#' @param n.obs the number of Gibbs samples.
#' @param treatment a set of given treatment assignment of length \code{m}.
#' @param covariates given confounder(s):
#' \itemize{
#' \item{\code{NULL}}{: no confounder.}
#' \item{a vector of length \code{m}}{: under unique confounder.}
#' \item{a \code{[q x m]} matrix}{: a set of \code{q} different confounders.}
#' }
#' 
#' @param initprob an initial probability generating outcomes. Defaults to \code{initprob} = 0.5
#' @param yvalues distinct binary values for outcomes. Defaults to \code{(0,1)}.
#' @param Neighborind a list of matrix specifying edge information of which first column illustrates a type of variables (1:outcome, 2:treatment, 3~:confounders) and of which second column presents the index of those variable.
#' @param Neighborpar index for parameters in the order of Neighborind.
#' @param n.burn the number of burn-in sample in Gibbs sampling (\eqn{\ge} \code{n.obs}). 
#' 
#' @return a \code{[n.obs x m]} matrix each row of which consists of outcomes. 
#' @export
#'
chaingibbs = function(pars, n.obs, treatment, covariates, initprob = 0.5, 
                      yvalues = c(0,1), Neighborind, Neighborpar,
                      n.burn){
  n.iter <- n.burn + n.obs
  m <- length(treatment)
 
  outcome <- matrix(0, n.iter, m)
  t <- 1
  outcome[t,] <- rbinom(m, 1, initprob)
  outcome[t,] <- ifelse(outcome[t,] == 1, max(yvalues), min(yvalues))
  observations <- rbind(outcome[t,], treatment, covariates)
  
  for (t in 2:n.iter) {
    random.index <- sample(1:m, 1)
    for (i in 1:m) {
      if (random.index == i) {
        out <- 0; out0 <- 0
        for (r in 1:length(Neighborpar[[i]])) {
          observations1 <- observations; observations1[1,i] <- max(yvalues)
          term <- as.matrix(observations1[Neighborind[[i]][[r]][,1], Neighborind[[i]][[r]][,2]])
          out <- out + pars[Neighborpar[[i]][[r]]]*prod(diag(term))
          observations0 <- observations; observations0[1,i] <- min(yvalues)
          term0 <- as.matrix(observations0[Neighborind[[i]][[r]][,1], Neighborind[[i]][[r]][,2]])
          out0 <- out0 + pars[Neighborpar[[i]][[r]]]*prod(diag(term0))
        }
      }
    }
    outcome[t,]  <- outcome[(t-1),]
    outcome[t, random.index] <- rbinom(1,1, exp(out)/(exp(out) + exp(out0)) )
    outcome[t, random.index] <- ifelse( outcome[t, random.index] == 1, max(yvalues), min(yvalues))
    observations[1,] <- outcome[t,]
  }
  
  Youtcomes <- outcome[c((n.burn+1):n.iter), ]
  return(Youtcomes)
}