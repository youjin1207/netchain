## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(netchain)
weight.matrix = matrix(c(0.5, 1, 0, 1, 0.3, 0.5, 0, 0.5, -0.5), 3, 3)
simobs = simGibbs(n.unit = 3, n.gibbs = 10, n.sample = 10, 
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
library(simcausal)
D <- DAG.empty()
D <- D +
  node("A1", t = 0, distr = "rbern", prob = 0.3) + 
  node("A2", t = 0, distr = "rbern", prob = 0.5) + 
  node("A3", t = 0, distr = "rbern", prob = 0.7) +
  node("Y1", t = 0, distr = "rbern", prob = plogis(-1 + A1[0])) + 
  node("Y2", t = 0, distr = "rbern", prob = plogis(0 - 0.5*A2[0])) + 
  node("Y3", t = 0, distr = "rbern", prob = plogis(1 - 0.5*A3[0]))
t.end <- 20
D <- D + 
  node("Y1", t = 1:t.end, distr = "rbern",
       prob = plogis(3*Y1[t-1] - 1.5 + A1[0] + Y2[t-1])) + 
  node("Y2", t = 1:t.end, distr = "rbern",
       prob = plogis(3*Y2[t-1] - 2 - 0.5*A2[0] + Y1[t-1] -Y3[t-1])) +
  node("Y3", t = 1:t.end, distr = "rbern", 
       prob = plogis(3*Y3[t-1] - 1 - 0.5*A3[0] -Y2[t-1]))

lDAG <- set.DAG(D)

## ------------------------------------------------------------------------
par(mar = c(4, 4, 1, 1))
plotDAG(lDAG, tmax = 10, xjitter = 0.3, yjitter = 0.1,
     edge_attrs = list(width = 0.5, arrow.width = 0.5, arrow.size = 0.5),
     vertex_attrs = list(size = 12, label.cex = 1, label.color = "dodgerblue"))

## ------------------------------------------------------------------------
par(mar = c(4, 4, 1, 1))
plotCG(lDAG)

