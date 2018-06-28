## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(netchain)
weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
simobs = simGibbs(n.unit = 3, n.gibbs = 100, n.sample = 10, 
                   weight.matrix,
                   treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
inputY = simobs$inputY
inputA = simobs$inputA
inputC = simobs$inputC
head(inputY)
head(inputA)
head(inputC)

## ------------------------------------------------------------------------
mat1 = rbind(c("Y", 1), c("A", 1), c("A", 2)) 

## ------------------------------------------------------------------------
mat2 = rbind(c("Y", 1), c("C1", 1), c("C2", 1)) 

## ------------------------------------------------------------------------
set.seed(1234)
library(netchain)
weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
simobs = simGibbs(n.unit = 3, n.gibbs = 100, n.sample = 10, 
                   weight.matrix,
                   treat.matrix = 0.5*diag(3), cov.matrix= (-0.3)*diag(3) )
inputY = simobs$inputY
inputA = simobs$inputA
inputC = simobs$inputC

R.matrix = ifelse(weight.matrix==0, 0, 1)      
diag(R.matrix) = 0

result = chain.causal.multi(targetoutcome = "mean", treatment = c(1,0,0), inputY, inputA, listC = inputC, R.matrix = R.matrix, E.matrix = diag(3), edgeinfo = list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), rbind(c("Y", 3), c("C", 3))), n.obs = 1000, n.burn = 100)

result

## ------------------------------------------------------------------------
set.seed(1234)
influence = causal.influence(targetoutcome = "mean", Avalues = c(1,0), 
                            inputY, inputA, listC = inputC, R.matrix, E.matrix = diag(3), 
                            edgeinfo = list(rbind(c("Y", 1), c("C", 1)), rbind(c("Y", 2), c("C", 2)), rbind(c("Y", 3), c("C", 3))), n.obs = 1000, n.burn = 100)
influence

