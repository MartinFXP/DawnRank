#' Dawn
#' 
#' The most basic version of DawnRank. This is DawnRank for one Patient. Dawn
#' differs from the other methods from calling DawnRank by being only for one
#' patient and not including the matrix information. Nonetheless it is still
#' quite useful for a quick hash for ranking a patient
#' 
#' 
#' @param dawnMatrix, the weighted adjacency
#' @param expressionVector, the normalized expression vector
#' @param mutationVector, a logical vector containing mutation information
#' @param damping, the damping vector
#' @param maxit, the maximum number of iterations to use, default 100
#' @param epsilon, the lower magnitude cutoff, default 0.0001
#' @param goldStandard, A list of common driver genes, used as a comparison.
#' This is optional, default=NULL
#' @param patientTag, an index for the patients
#' @return the ranks. A list of 3 including a [[1]] output of all the ranks,
#' [[2]] mutated ranks, [[3]] the steps of convergence
#' @export Dawn
Dawn <-
function(dawnMatrix, expressionVector, mutationVector,damping, maxit=100, epsilon=0.0001, goldStandard=NULL,patientTag="defaultPatient"){
  
  ## dawnMatrix<-DawnMatrix(adjMatrix)
  
  #need to make a column matrix for the ranking t
  ranking_T<-matrix(expressionVector/sum(expressionVector),length(expressionVector),1)
  #  ranking_T<-matrix(expressionVector,length(expressionVector),1)
  
  #there is a constant term and a variable term in the equation. The constant term is ex*damping, and the iterative term is everything else
  constantTerm<-expressionVector/sum(expressionVector)*(1-damping)
  #  constantTerm<-expressionVector*(1-damping)
  
  #do the iterations until convergence
  for(i in 1:maxit){
    #matrix to get the summation of the rank, this is just a simplification of the summation of the ranks
    ranking_T1<-damping*dawnMatrix%*%ranking_T+constantTerm        
    mag<-sqrt(sum((ranking_T1-ranking_T)^2))
    #update the Rank table
    ranking_T<-ranking_T1
    if(mag < epsilon) break
  }
  message(paste("completed at iteration i:",i))
  
  #this is to organize the data
  outFrame<-data.frame(Rank=ranking_T,PercentRank=100*rank(ranking_T)/(length(ranking_T)+1),isMutated=mutationVector)
  if(length(goldStandard)>0){
    outFrame$isGoldStandard<-0
    outFrame$isGoldStandard[match(goldStandard,rownames(outFrame))]<-1
  }
  #find only the summary of mutated genes
  if(sum(outFrame$isMutated)>0){
    mutatedRanks<-outFrame[outFrame$isMutated==1,]
    mutatedRanks<-mutatedRanks[,-3]
    mutatedRanks<-data.frame(Gene=rownames(mutatedRanks),Patient=patientTag,mutatedRanks)
    rownames(mutatedRanks)<-NULL
  }else{
    mutatedRanks<-logical(0)
  }
  
  return(list(summaryOutput=outFrame,mutatedRanks=mutatedRanks,convergenceIterations=i))
  
}
