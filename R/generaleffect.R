#' Causal estimation on collective outcomes under multiple confounders and interference.
#' 
#' This function calculates probability associated with counterfactual collective outcome(s)
#'  P(\strong{Y}(\strong{a}) = \strong{y}) when \code{m} units are subject to interference 
#' and contagion possibly with the presence of multiple confounders. To estimate the magnitude of 
#' main effects, two-way interaction effects, or any higher-order interaction effects we use hybrid
#' graphcial models combining features of both log-linear models on undirected graphs (\code{R.matrix}) 
#' and directed acyclic graphs (DAGs) models used to represent casual relationships. 
#' 
#' 
#'
#' @param targetoutcome is a targeted couterfactual outcome of probability is in our interest, having different forms: 
#' \describe{
#'    \item{a vector of length \code{m}}{a vector specifies every element of \strong{y}.}
#'    \item{a \code{[q x m]} matrix}{a collection of \strong{y_1}, \strong{y_2}, ..., \strong{y_q} of which we want to derive the probability.}
#'    \item{an integer}{the number of 1's in \strong{y} (\eqn{0 \ge} & \eqn{\le m}).}
#'    \item{'mean'}{when we want derive E(\strong{Y}(\strong{a})) (default).}
#' }
#' @param treatment a vector of length \code{m} representing given treatment assignment \strong{a}.
#' @param inputY a \code{[n x m]} matrix of \code{n} independent outcomes for \code{m} units.
#' @param inputA a \code{[n x m]} matrix of \code{n} independent treatment assignments assigned to \code{m} units.
#' @param listC is either a matrix, list or \code{NULL}:
#' \describe{
#'    \item{a \code{[n x m]} matrix}{a matrix of \code{n} independent confounders for \code{m} units under single confounder.}
#'    \item{a list of \code{[n x m]} matrices}{a collection of \code{n} independent confounders for \code{m} units under multiple confounders.}
#'    \item{\code{NULL}}{no confounders.}
#' }
#' @param R.matrix a \code{[m x m]} relational symmetric matrix where \eqn{R.matrix_ij = 1} indicates \eqn{Y_i} and \eqn{Y_j} are adjacent. 
#' @param E.matrix a \code{[m x m]} matrix where \eqn{E.matrix_ij = 1} indicates \eqn{A_i} has a direct causal effect on \eqn{Y_j}. Defaults to diagonal matrix, which indicates no interference. 
#' @param edgeinfo a list of matrix specifying additional directed edges (from confounders or treatment to the outcomes) information. Defaults to \code{NULL}.
#' \describe{
#'    \item{first column:}{\code{"Y"} indicates outcomes; \code{"A"} indicates treatment; \code{"C"} indicates confounders. Under multiple confounders, \code{"C1"}, \code{"C2"}, ... indicate each confounder.}
#'    \item{second column:}{an index for unit corresponding to the variable in the first column, \code{i=1,2,...m}.}
#' }
#' @param n.obs the number of Gibbs samplers except for burn-in sample.
#' @param n.burn the number of burn-in sample in Gibbs sampling.
#' @param optim.method the method used in \code{optim()}. Defaults to \code{"L-BFGS-B"}.
#'
#' @return returns \code{"noconvergence"} in case of failure to converence or a list with components :
#' \item{\code{causalprob}}{the estimated probability P(\strong{Y}(\strong{a}) = \strong{y}).}
#' \item{\code{n.par}}{the number of parameters estimated in conditional log-linear model.}
#' \item{\code{par.est}}{the estimated parameters.}
#' 
#' @author Youjin Lee
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom gtools permutations
#' @importFrom stringr str_extract
#' @importFrom stats optim 
#' @importFrom stats rbinom
#' 
#' @export 
#'
#' @examples
#' library(netchain)
#' set.seed(1234)
#' weight.matrix <- matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
#' simobs <- simGibbs(n.unit = 3, n.gibbs = 100, n.sample = 5, 
#'                   weight.matrix, treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
#' inputY <- simobs$inputY                   
#' inputA <- simobs$inputA   
#' inputC <- simobs$inputC 
#' R.matrix <- ifelse(weight.matrix==0, 0, 1)    
#' diag(R.matrix) <- 0
#' edgeinfo <- list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), 
#'            rbind(c("Y", 3), c("C", 3)))  
#' # implement a function (take > 10 seconds)
#' # result <- chain.causal.multi(targetoutcome = "mean",
#' # treatment <- c(1,0,0), inputY, inputA, listC = inputC, R.matrix, 
#' # E.matrix <- diag(3), edgeinfo = edgeinfo)
#' 
#' 
chain.causal.multi = function(targetoutcome = "mean", treatment, inputY, inputA, listC, R.matrix, E.matrix, edgeinfo = NULL, 
                              n.obs = 1000, n.burn = 100, optim.method = "L-BFGS-B"){
  
  allobservations = list()
  for(i in 1:nrow(inputY)){
    allobservations[[i]] = rbind(inputY[i,], inputA[i,])
    if(class(listC) == "list"){
      for(r in 1:length(listC)){
        allobservations[[i]] = rbind(allobservations[[i]], listC[[r]][i,])
      }
    }else if(class(listC) == "matrix"){
      allobservations[[i]] = rbind(allobservations[[i]], listC[i,])
    }
  }
  
  yvalues = unique(as.numeric(inputY))
  
  if(!is.null(edgeinfo)){
    edgeExtra = list()
    for(k in 1:length(edgeinfo)){
      tmpname = ifelse(edgeinfo[[k]][,1] == "Y", 1, 0)
      tmpname = ifelse(edgeinfo[[k]][,1] == "A", 2, tmpname)
      if(sum(str_extract(edgeinfo[[k]][,1], "[aA-zZ]+") == "C")!=0){
        var.order =  which(str_extract(edgeinfo[[k]][,1], "[aA-zZ]+") == "C")
        confounder.num = str_extract(edgeinfo[[k]][var.order,1], "[0-9]+")
        if(is.na(confounder.num)) confounder.num = 1
        tmpname[var.order] = 2 + as.integer(confounder.num)
      }
      edgeExtra[[k]] = cbind(as.integer(tmpname), as.integer(edgeinfo[[k]][,2]))
    }
  }else{
    edgeExtra = NULL
  }
  
  # define edgeY; edgeAY; edgeCY
  edgeY = cbind(which(R.matrix == 1) %/% nrow(R.matrix) + 1, which(R.matrix == 1) %% nrow(R.matrix))
  edgeY[which(edgeY[,2] == 0), 1] = edgeY[which(edgeY[,2] == 0),1] -1
  edgeY[which(edgeY[,2] == 0), 2] = nrow(R.matrix)
  edgeY = edgeY[which(edgeY[,1] < edgeY[,2]),]
  if(class(edgeY) == "numeric") edgeY = t(as.matrix(edgeY))
  colnames(edgeY) = c("Y", "Y")
  
  edgeAY = cbind(which(E.matrix == 1) %% nrow(E.matrix), which(E.matrix == 1) %/% nrow(E.matrix) + 1)
  edgeAY[which(edgeAY[,1] == 0), 2] = edgeAY[which(edgeAY[,1] == 0),2] -1
  edgeAY[which(edgeAY[,1] == 0), 1] = nrow(E.matrix)
  if(class(edgeAY) == "numeric") edgeAY = t(as.matrix(edgeAY))
  colnames(edgeAY) = c("A", "Y")
  
  ## define n.par and par.est ##
  n.par = ncol(inputY) + nrow(edgeY) + nrow(edgeAY) + length(edgeExtra)
  permutetab = permutations(n=length(unique(as.numeric(inputY))), r=ncol(inputY), unique(as.numeric(inputY)), repeats.allowed=T)
  
  par.est = try(optim(par = rep(0, n.par), multiloglikechain, listobservations = allobservations,
                      permutetab = permutetab, edgeY = edgeY,
                      edgeAY = edgeAY, edgeExtra = edgeExtra,
                      control = list(fnscale = -1), method = optim.method)$par,
                silent = TRUE)
  if(class(par.est) == "try-error") return("noconvergence")
  
  #multiloglikechain(NumericVector pars, List allobservations, NumericMatrix permutetab)
  Neighborind = Neighborpar = list()
  for(i in 1:ncol(inputY)){
    Neighborind[[i]] = list(); Neighborpar[[i]] = list();
    Neighborind[[i]][[1]] = t(as.matrix(c(1,i))); Neighborpar[[i]][[1]] = i;
    whichnb = which(rowSums(edgeY == i)!=0)
    if(length(whichnb) > 0){
      for(j in 1:length(whichnb)){
        Neighborind[[i]][[1+j]] = rbind(c(1,i) ,cbind(1, edgeY[whichnb[j],][edgeY[whichnb[j],]!=i])); 
        Neighborpar[[i]][[1+j]] = ncol(inputY) + whichnb[j];
      }
    }
    whicheffect = which(edgeAY[,2] == i) # first column : A; second column : Y
    if(length(whicheffect) > 0){
      for(l in 1:length(whicheffect)){
        Neighborind[[i]][[1+length(whichnb)+l]] = rbind(c(1,i) ,cbind(2, edgeAY[whicheffect[l],1])); 
        Neighborpar[[i]][[1+length(whichnb)+l]] = ncol(inputY) + nrow(edgeY) + whicheffect[l];
      }
    }
    count = 0
    for(k in 1:length(edgeExtra)){
      if(sum(edgeExtra[[k]][,1] == 1 & edgeExtra[[k]][,2] == i) > 0){
        count = count + 1
        mynode = which(edgeExtra[[k]][,1] == 1 & edgeExtra[[k]][,2] == i)
        Neighborind[[i]][[1+length(whichnb)+length(whicheffect)+count]] = rbind(c(1,i) ,cbind(edgeExtra[[k]][-mynode,1], edgeExtra[[k]][-mynode,2]))
        Neighborpar[[i]][[1+length(whichnb)+length(whicheffect)+count]] = ncol(inputY) + nrow(edgeY) + nrow(edgeAY) + k
      }
    }
  }
  
  
  ## g-formula under the existence of confounders
  targets = 0
  for(i in 1:length(allobservations)){
    if(nrow(allobservations[[1]]) > 2){
      covariates = allobservations[[i]][3:nrow(allobservations[[i]]),]
    }else{
      covariates = NULL
    }
    outcomes =  chaingibbs(pars = par.est, n.obs = n.obs, treatment, covariates, initprob = 0.5, yvalues = c(0,1), Neighborind, Neighborpar,
                           n.burn = n.burn)
    if(class(targetoutcome) == "numeric" & length(targetoutcome) == ncol(inputY)){
      targets = targets + mean(rowMeans(outcomes == targetoutcome) == 1 ) / length(allobservations)
    }else if(class(targetoutcome) == "matrix" ){
      for(jj in 1:nrow(targetoutcome)){
        targets = targets +  mean(apply(outcomes, 1, function(x) identical(x, targetoutcome[jj,])))  / length(allobservations)
      }
    }else if(class(targetoutcome) == "numeric" & length(targetoutcome) == 1){
      targets = targets +  mean(rowSums(outcomes == max(yvalues)) == targetoutcome ) / length(allobservations)
    }else{
      targets = targets + mean(rowMeans(outcomes == max(yvalues))) / length(allobservations)
    }
  }
  return(list(causalprob = targets, n.par = n.par, par.est = par.est))
}