#' Generate binary (\strong{Y}, \strong{A}, \strong{C}) from chain graph model under simplest scenario. 
#'
#' @param n.unit the number of units (\code{m}).
#' @param n.obs the number of independent observations (\code{n}).
#' @param weight.matrix a \code{[m x m]} symmetric relational matrix where \eqn{W_ij = 1} indicates the existence of undirected edges between \eqn{Y_i} and \eqn{Y_j} and its magnitude. Here \eqn{W_ii} represents the main effect of \eqn{Y_i}.  
#' @param treat.matrix a \code{[m x m]} matrix where \eqn{treat.matrix_ij} indicates the magnitude of direct effect from \eqn{A_i} to \eqn{Y_j}.
#' @param cov.matrix a \code{[m x m]} matrix where \eqn{treat.matrix_ij} indicates the magnitude of direct effect from \eqn{C_i} to \eqn{Y_j}.
#' @param init.prob an initial probability generating (\strong{Y}, \strong{A}, \strong{C}) from Bernoulli distribution.
#' @param treat.prob a probability updating \strong{A} in Gibbs sampling procedure.
#' @param cov.prob a probability updating \strong{C} in Gibbs sampling procedure. 
#' @param n.burn the number of burn-in sample in Gibbs sampling (\eqn{\ge} \code{n.obs}). 
#'
#' @return a list of \code{[n.obs]} independent observations:
#' \item{\code{inputY}}{a \code{[n.obs x m]} matrix for outcomes.}
#' \item{\code{inputA}}{a \code{[n.obs x m]} matrix for treatments.}
#' \item{\code{inputC}}{a \code{[n.obs x m]} matrix for confounders.}
#' 
#' @export
#'
#' @examples
#' library(netchain)
#' weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
#' simobs = simGibbs(n.unit = 3, n.obs = 200, 
#'                   weight.matrix,
#'                   treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
#' result = simobs
#' inputY = result$inputY
#' inputA = result$inputA
#' inputC = result$inputC
#' 
simGibbs = function(n.unit, n.obs, 
                    weight.matrix, treat.matrix, cov.matrix,
                    init.prob = 0.5, treat.prob = 0.5, cov.prob = 0.5, 
                    n.burn = 100){
  
  n.iter = n.burn + n.obs
  treatment = cov = outcome = matrix(0, n.iter, n.unit)
  t = 1
  treatment[t,] = rbinom(n.unit, 1, init.prob)
  cov[t,] = rbinom(n.unit, 1, init.prob)
  outcome[t,]  = rbinom(n.unit, 1, init.prob)
  
  for(t in 2:n.iter){
    random.ind = sample(1:n.unit, 1)
    # from previous time
    treatment[t,] = treatment[t-1, ] 
    cov[t,] = cov[t-1,]
    outcome[t,] = outcome[t-1,]
    # update new treatment and covariate
    cov[t,random.ind] = rbinom(1,1,cov.prob)
    treatment[t,random.ind] = rbinom(1,1,treat.prob)
    
    for(i in 1:n.unit){
      if(random.ind == i){
        out =  sum(outcome[t-1,]*weight.matrix[i,]) +
          sum(treatment[t-1,]*treat.matrix[,i]) +
          sum(cov[t-1,]*cov.matrix[i,]) 
      }
    }
    outcome[t, random.ind] = rbinom(1,1, exp(out) / (1 + exp(out)))
  }
  
  inputY = outcome[c((n.burn+1):n.iter), ]
  inputA = treatment[c((n.burn+1):n.iter), ]
  inputC = cov[c((n.burn+1):n.iter), ]
  
  return(list(inputY = inputY, inputA = inputA, inputC = inputC))
  
}