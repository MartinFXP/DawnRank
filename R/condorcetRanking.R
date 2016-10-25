#' condorcetRanking
#' 
#' condorcetRanking determines the aggregate rank of a gene in when given
#' population data
#' 
#' 
#' @param scoreMatrix, a matrix containing all the given DawnRank scores per
#' patient. Rows are genes, columns are patients
#' @param mutationMatrix, the mutation matrix. Rows are genes, columns are
#' patients
#' @param pen, the penality parameter in the condorcet algorithm for missing
#' data. Default 0.85
#' @return the ranks. A list of 2 including a [[1]] a matrix of all pairwise
#' comparisons, [[2]] the final rankings based on the Condorcet score
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
#' @export condorcetRanking
condorcetRanking <-
function(scoreMatrix,mutationMatrix,pen=0.85){
  pretrunc<-apply(mutationMatrix,1,sum)
  trunc<-names(pretrunc[pretrunc>0])
  trunc<-intersect(trunc,rownames(scoreMatrix))
  pwc<-function(p1,p2,scoreMatrix,mutationMatrix){
    voters<-colnames(mutationMatrix)[(mutationMatrix[p1,]!=0)|(mutationMatrix[p2,]!=0)]
    fightmat<- as.matrix(scoreMatrix[c(p1,p2),voters])
    p1wins<-sum((fightmat[1,]-fightmat[2,])>0)
    p2wins<-ncol(fightmat)-p1wins
    return(c(p1wins,p2wins))
  }
  scoreMatrix[mutationMatrix==0]<-pen*scoreMatrix[mutationMatrix==0]
  scoreMatrix<-scoreMatrix[trunc,]
  cmat<-matrix(0,nrow=length(trunc),ncol=length(trunc),dimnames=list(trunc, trunc))
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      if(j>i){
        vres<-pwc(trunc[i],trunc[j],scoreMatrix,mutationMatrix)
        cmat[i,j]<-vres[1]
        cmat[j,i]<-vres[2]
      }
    }
  }
  
  wins<-apply(cmat,1,sum)
  losses<-apply(cmat,2,sum)
  copelandcriterion<-sort(wins/(wins+losses),decreasing=TRUE)
  return(list(cmat,copelandcriterion))
}
