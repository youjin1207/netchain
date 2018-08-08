#' Generate binary (\strong{Y}, \strong{A}, \strong{C}) from chain graph model under simplest scenario. 
#'
#' @param n.unit the number of units (\code{m}).
#' @param n.gibbs the number of independent Gibbs sampler.
#' @param n.sample the number of samples from each Gibbs sampling (\code{n} = \code{n.gibbs} x \code{n.sample}). 
#' @param weight.matrix a \code{[m x m]} symmetric relational matrix where \eqn{W_ij = 1} indicates the existence of undirected edges between \eqn{Y_i} and \eqn{Y_j} and its magnitude. Here \eqn{W_ii} represents the main effect of \eqn{Y_i}.  
#' @param treat.matrix a \code{[m x m]} matrix where \eqn{treat.matrix_ij} indicates the magnitude of direct effect from \eqn{A_i} to \eqn{Y_j}.
#' @param cov.matrix a \code{[m x m]} matrix where \eqn{treat.matrix_ij} indicates the magnitude of direct effect from \eqn{C_i} to \eqn{Y_j}.
#' @param init.prob an initial probability generating (\strong{Y}, \strong{A}, \strong{C}) from Bernoulli distribution.
#' @param treat.prob a probability updating \strong{A} in Gibbs sampling procedure.
#' @param cov.prob a probability updating \strong{C} in Gibbs sampling procedure. 
#' @param n.burn the number of burn-in sample in Gibbs sampling (\eqn{\ge} \code{n.obs}). 
#' @param yvalues a vector of distinct binary outcome values. Defaults to \code{c(0,1)}.
#'
#' @return a list of \code{[n.gibbs] x [n.sample]} independent observations:
#' \item{\code{inputY}}{a \code{[([n.gibbs] x [n.sample]) x m]} matrix for outcomes.}
#' \item{\code{inputA}}{a \code{[([n.gibbs] x [n.sample]) x m]} matrix for treatments.}
#' \item{\code{inputC}}{a \code{[([n.gibbs] x [n.sample]) x m]} matrix for confounders.}
#' 
#' @export
#'
#' @examples
#' library(netchain)
#' weight.matrix <- matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
#' simobs <- simGibbs(n.unit = 3, n.gibbs = 200, n.sample = 10, 
#'                   weight.matrix,
#'                   treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
#' inputY <- simobs$inputY
#' inputA <- simobs$inputA
#' inputC <- simobs$inputC
#' 
simGibbs <- function(n.unit, n.gibbs, n.sample, 
                    weight.matrix, treat.matrix, cov.matrix,
                    init.prob = 0.5, treat.prob = 0.5, cov.prob = 0.5, 
                    n.burn = 100, yvalues = c(1,0)){
  inputY = inputA = inputC = c()
  
  for (r in 1:n.gibbs) {
    n.iter <- n.burn + n.sample
    outcome <- matrix(0, n.iter, n.unit)
    t <- 1
    treatment <- rbinom(n.unit, 1, init.prob)
    cov <- rbinom(n.unit, 1, init.prob)
    outcome[t,]  <- rbinom(n.unit, 1, init.prob)
    outcome[t,] <- ifelse(outcome[t,] == 1, max(yvalues), min(yvalues))
    for (t in 2:n.iter) {
      random.ind <- sample(1:n.unit, 1)
      # from previous time
      outcome[t,] <- outcome[t-1,]
      for (i in 1:n.unit) {
        if (random.ind == i) {
          out <- diag(weight.matrix)[i] +  sum(outcome[t-1,-i]*weight.matrix[i,-i]) +
            sum(treatment*treat.matrix[,i]) + sum(cov*cov.matrix[i,]) 
        }
      }
      outcome[t, random.ind] <- rbinom(1,1, exp(out) / (1 + exp(out)))
      outcome[t, random.ind] <- ifelse(outcome[t, random.ind] == 1, max(yvalues), min(yvalues))
    }
    
    inputY <- rbind(inputY, outcome[c((n.burn+1):n.iter),])
    inputA <- rbind(inputA, t(matrix(rep(treatment, n.sample), n.unit, n.sample)))
    inputC <- rbind(inputC, t(matrix(rep(cov, n.sample), n.unit, n.sample)))
  }
  return(list(inputY = inputY, inputA = inputA, inputC = inputC))
}