#' Plot chain graph components from DAG
#'
#' @param DAG A DAG object from \code{simcausal} R Package.
#' @param vertex_attrs A list of parameters for DAG vertices which are passed on to \code{add.vertices} in \code{igraph}.
#' @param edge_attrs  A list of parameters for DAG edges which are passed on to \code{add.edges} in \code{igraph}.
#'
#' @return plot
#' @export
#' 
#' @import simcausal
#' @import igraph
#'
#' @examples
#' library(netchain)
#' library(simcausal)
#' D <- DAG.empty()
#' D <- D +
#' node("A1", t = 0, distr = "rbern", prob = 0.3) + 
#' node("A2", t = 0, distr = "rbern", prob = 0.5) + 
#' node("A3", t = 0, distr = "rbern", prob = 0.7) +
#' node("Y1", t = 0, distr = "rbern", prob = plogis(-1 + A1[0])) + 
#' node("Y2", t = 0, distr = "rbern", prob = plogis(0 - 0.5*A2[0])) + 
#' node("Y3", t = 0, distr = "rbern", prob = plogis(1 - 0.5*A3[0])) + 
#' node("Y4", t = 0, distr = "rbern", prob = plogis(1 - 0.5*A3[0]))
#' t.end <- 20
#' D <- D + 
#' node("Y1", t = 1:t.end, distr = "rbern",
#' prob = plogis(3*Y1[t-1] - 1.5 + A1[0] + Y2[t-1])) + 
#' node("Y2", t = 1:t.end, distr = "rbern",
#' prob = plogis(3*Y2[t-1] - 2 - 0.5*A1[0] + Y1[t-1])) +
#' node("Y3", t = 1:t.end, distr = "rbern", 
#' prob = plogis(3*Y3[t-1] - 1 - 0.5*A2[0] -Y4[t-1])) + 
#' node("Y4", t = 1:t.end, distr = "rbern", 
#' prob = plogis(3*Y4[t-1] - 1 - 0.5*A3[0] -Y3[t-1]))
#' lDAG <- set.DAG(D)
#' plotCG(lDAG, vertex_attrs = list(label.color = "hotpink"))
#' 
#' 
plotCG <- function(DAG, vertex_attrs = list(), edge_attrs = list()){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph R package is required.", call. = FALSE)
  }
  
  ### variable without no parents
  parents <- c()
  for(i in 1:length(names(DAG))){
    if (length(attr(DAG, "parent")[[i]]) == 0) parents <- c(parents, names(DAG)[[i]])
  }
  
  max.time <- max(as.integer(unlist(lapply( strsplit(names(DAG), "_"), '[[', 2))))
  
  varnames <- unlist(lapply( strsplit(names(DAG), "_"), '[[', 1)) 
  time <- unlist(lapply( strsplit(names(DAG), "_"), '[[', 2)) 
  chainset <- unique(varnames[as.integer(time) >= 2])
  if (length(chainset) == 0) return("No chain component detected.")
  parentset = friendset = list()
  
  for (i in 1:length(chainset)) {
    
    t <- 1
    tmp <- eval(parse(text = as.character(paste0("attr(DAG, \"parent\")", "$", chainset[i], "_" ,as.character(t)))))
    parentset[[i]] <- parents[parents %in% tmp]
    friendset[[i]] <-  unlist(lapply(strsplit(tmp, "_"), '[[', 1 ))[which(as.integer(unlist(lapply(strsplit(tmp, "_"), '[[', 2 ))) == t-1)]
    
    for (t in 2:max.time) {
      tmp <- eval(parse(text = as.character(paste0("attr(DAG, \"parent\")", "$", chainset[i], "_" ,as.character(t)))))
      
      newparentset <- parents[parents %in% tmp] 
      parentset[[i]] <- newparentset[newparentset %in% parentset[[i]]]
      
      newfriendset <- unlist(lapply(strsplit(tmp, "_"), '[[', 1 ))[which(as.integer(unlist(lapply(strsplit(tmp, "_"), '[[', 2 ))) == t-1)]
      friendset[[i]] <- newfriendset[newfriendset %in% friendset[[i]]]
    }
  }
  parentnode <- unique(unlist(parentset))
  friendnode <- unique(unlist(friendset))
  
  edgelist <- c()
  for (i in 1:length(chainset)) {
    for (r in 1:length(parentset[[i]])) {
      edgelist <- rbind(edgelist, c(chainset[i], parentset[[i]][[r]]))
    }
    for (r in 1:length(friendset[[i]])) {
      edgelist <- rbind(edgelist, c(chainset[i], friendset[[i]][[r]]))
    }
  }
  ## check validity of chain graph
  deleterow <- c()
  for (i in 1:length(chainset)) {
    if (sum(rowSums(edgelist == chainset[i]) == 2) == 0) return("Invalid chain graph approximation.")
    deleterow <- c(deleterow, which(rowSums(edgelist == chainset[i]) == 2))
  }
  if (!is.null(deleterow)) edgelist = edgelist[-deleterow,]
  
  deleterow <- c()
  for (i in 1:nrow(edgelist)) {
    if (sum(edgelist[i,] %in% chainset) == 2) {
      if (sum(rowSums(matrix(edgelist %in% edgelist[i,], ncol = 2)) == 2 ) < 2) deleterow <- c(deleterow, i) 
    }
  }
  if (!is.null(deleterow)) edgelist <- edgelist[-deleterow,]
  
  ## draw graph
  attnames_ver <- names(vertex_attrs)
  attnames_edge <- names(edge_attrs)
  vertex_attrs_default <- list(color = NA, label.color = "seagreen", shape = "none", size = 15, label.cex = 1, label.dist = 0)
  edge_attrs_default <- list(color = "black", width = 1, lty = 1, arrow.width = 1, arrow.size = 1)
  vertex_attrs <- append(vertex_attrs, vertex_attrs_default[!(names(vertex_attrs_default)%in%attnames_ver)])
  edge_attrs <- append(edge_attrs, edge_attrs_default[!(names(edge_attrs_default)%in%attnames_edge)])
  ####
  
  g <- igraph::graph.empty()
  g <- igraph::add.vertices(g, nv = length(parentnode) + length(friendnode), attr = vertex_attrs)
  igraph::V(g)$name <- c(parentnode, friendnode)
  g <- igraph::add.edges(g, as.vector(t(edgelist)), attr = edge_attrs)
  
  
  #g <- igraph::set.graph.attribute(g, 'layout')
  edgemode <- ifelse(rowSums(matrix(as_edgelist(g) %in% parentnode, ncol = 2) > 0), 1, 0)
  igraph::plot.igraph(g, edge.arrow.mode = edgemode)
  
}
