\name{DawnRank}
\alias{DawnRank}
\title{DawnRank}
\usage{
  DawnRank(adjMatrix, expressionMatrix, mutationMatrix,
    mu = 20, maxit = 100, epsilon = 1e-04,
    goldStandard = NULL)
}
\arguments{
  \item{adjMatrix,}{the adjacency matrix}

  \item{expressionMatrix,}{the normalized expression matrix
  (multiple patients)}

  \item{mutationMatrix,}{a logical matrix containing
  mutation information}

  \item{mu,}{the proposed free parameter}

  \item{maxit,}{the maximum number of iterations to use,
  default 100}

  \item{epsilon,}{the lower magnitude cutoff, default
  0.0001}

  \item{goldStandard,}{A list of common driver genes, used
  as a comparison. This is optional, default=NULL}
}
\value{
  the ranks. A list of 3 including a [[1]] output of all
  the ranks, [[2]] output of all the ranks (percentile),
  [[3]] mutated ranks, [[4]] the steps of convergence
}
\description{
  This is the main method of the DawnRank method. DawnRank calculates the importance of a gene in a pathway through its connectivity to downstream genes and the differential expression of the downstream genes.  
}


\examples{
###using a small subset of the TCGA dataset and a small KEGG 
###gene interaction network,
###We will show how to get DawnRank Results

library(DawnRank)

#load the mutation data
data(brcaExampleMutation)

#load the tumor expression data
data(brcaExampleTumorExpression)

#load the normal expression data
data(brcaExampleNormalExpression)

#load the pathway data
data(brcaExamplePathway)

#load the gold standard
data(goldStandard)

#normalize the tumor and normal data to get the differential expression
normalizedDawn<-DawnNormalize(tumorMat=brcaExampleTumorExpression,
normalMat=brcaExampleNormalExpression)

#get the DawnRank Score
dawnRankScore<-DawnRank(adjMatrix=brcaExamplePathway,
mutationMatrix=brcaExampleMutation,expressionMatrix=normalizedDawn, 
mu=3,goldStandard=goldStandard)

#look at the DawnRank scores for a few patients
dawnRankFrame<-dawnRankScore[[3]]
head(dawnRankFrame)

#get the aggregate DawnRank scores
aggregateDawnRankScore<-condorcetRanking(scoreMatrix=dawnRankScore[[2]],
mutationMatrix=brcaExampleMutation)

#look at top 10 ranked genes
top10<-aggregateDawnRankScore[[2]][1:10]
top10

#get the individual cutoff for patient TCGA-A2-A04P
dawnRankFrame$isCGC<-dawnRankFrame$isGoldStandard
library(maxstat)
significance<-patspeccutoff(patient="TCGA-A2-A04P",ms=dawnRankFrame,
default=95)
}
