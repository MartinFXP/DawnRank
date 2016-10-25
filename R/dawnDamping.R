#' dawnDamping
#' 
#' dawnDamping calculates the damping factor
#' 
#' 
#' @param adjMatrix, the adjacency matrix
#' @param mu, the free parameter to fit
#' @return the vector for individualized damping factors
#' @export dawnDamping
dawnDamping <-
function(adjMatrix,mu){
  links<-apply(adjMatrix,2,sum)
  dvec<-links/(links+mu)
  return(dvec)
}
