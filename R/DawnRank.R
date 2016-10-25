#' DawnRank
#' 
#' DawnRank is a general version that can analyze many patients at once
#' 
#' 
#' @param adjMatrix, the adjacency matrix
#' @param expressionMatrix, the normalized expression matrix (multiple
#' patients)
#' @param mutationMatrix, a logical matrix containing mutation information
#' @param mu, the proposed free parameter
#' @param maxit, the maximum number of iterations to use, default 100
#' @param epsilon, the lower magnitude cutoff, default 0.0001
#' @param goldStandard, A list of common driver genes, used as a comparison.
#' This is optional, default=NULL
#' @return the ranks. A list of 3 including a [[1]] output of all the ranks,
#' [[2]] output of all the ranks (percentile), [[3]] mutated ranks, [[4]] the
#' steps of convergence
#' @examples
#' 
#' ###using a small subset of the TCGA dataset and a small KEGG gene interaction network,
#' ###We will show how to get DawnRank Results
#' 
#' library(DawnRank)
#' 
#' #load the mutation data
#' data(brcaExampleMutation)
#' 
#' #load the tumor expression data
#' data(brcaExampleTumorExpression)
#' 
#' #load the normal expression data
#' data(brcaExampleNormalExpression)
#' 
#' #load the pathway data
#' data(brcaExamplePathway)
#' 
#' #load the gold standard
#' data(goldStandard)
#' 
#' #normalize the tumor and normal data to get the differential expression
#' normalizedDawn<-DawnNormalize(tumorMat=brcaExampleTumorExpression,normalMat=brcaExampleNormalExpression)
#' 
#' #get the DawnRank Score
#' dawnRankScore<-DawnRank(adjMatrix=brcaExamplePathway,mutationMatrix=brcaExampleMutation,expressionMatrix=normalizedDawn, mu=3,goldStandard=goldStandard)
#' 
#' #look at the DawnRank scores for a few patients
#' dawnRankFrame<-dawnRankScore[[3]]
#' head(dawnRankFrame)
#' 
#' #get the aggregate DawnRank scores
#' aggregateDawnRankScore<-condorcetRanking(scoreMatrix=dawnRankScore[[2]],mutationMatrix=brcaExampleMutation)
#' 
#' #look at top 10 ranked genes
#' top10<-aggregateDawnRankScore[[2]][1:10]
#' top10
#' 
#' #get the individual cutoff for patient TCGA-A2-A04P
#' dawnRankFrame$isCGC<-dawnRankFrame$isGoldStandard
#' library(maxstat)
#' patspeccutoff(patient="TCGA-A2-A04P",ms=dawnRankFrame,default=95)->significance
#' 
#' @export DawnRank
DawnRank <-
function(adjMatrix, expressionMatrix, mutationMatrix, mu=20, maxit=100, epsilon=0.0001, goldStandard=NULL){
  damping<-dawnDamping(adjMatrix,mu)
  allRanks<-logical(0);allPercentiles<-logical(0);allMutated<-logical(0); iterations<-logical(0)
  #run Dawn through a loop
  for(i in 1:ncol(expressionMatrix)){
    Dawni<-Dawn(adjMatrix,expressionMatrix[,i],mutationMatrix[,i],damping=damping,maxit=maxit,epsilon=epsilon,patientTag=colnames(expressionMatrix)[i],goldStandard=goldStandard) 
    allPercentiles<-cbind(allPercentiles,Dawni[[1]][,2])
    allRanks<-cbind(allRanks,Dawni[[1]][,1])
    allMutated<-rbind(allMutated,Dawni[[2]])
    iterations<-append(iterations,Dawni[[3]])
  }
  rownames(allPercentiles)<-rownames(expressionMatrix)
  colnames(allPercentiles)<-colnames(expressionMatrix)
  rownames(allRanks)<-rownames(expressionMatrix)
  colnames(allRanks)<-colnames(expressionMatrix)
  allMutated$Frequency<-table(allMutated[,1])[allMutated[,1]]
  averageRank<-apply(allRanks,1,mean)
  sdRank<-apply(allRanks,1,sd)
  allMutated$Deviation<-(allMutated$Rank-averageRank[paste(allMutated$Gene)])/sdRank[paste(allMutated$Gene)]
  
  return(list(AllRanks=allRanks,AllPercentiles=allPercentiles,AllMutated=allMutated,iterations=iterations))
}
