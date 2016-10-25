#' DawnMatrix
#' 
#' DawnMatrix creates a transition matrix for DawnRank
#' 
#' 
#' @param adjMatrix, the adjacency matrix
#' @return the transition matrix
#' @export DawnMatrix
DawnMatrix <-
function(adjMatrix){
  
  #get the sums of the columns of the adjacency matrix
  colsums<-apply(adjMatrix,2,sum)
  transitionMatrix<-t(t(adjMatrix)/(colsums+1e-16))
  return(transitionMatrix)
}
