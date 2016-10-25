\name{DawnNormalize}
\alias{DawnNormalize}
\title{DawnNormalize}
\usage{
  DawnNormalize(tumorMat, normalMat)
}
\arguments{
  \item{tumorMat,}{a matrix representing the tumor
  expression}

  \item{normalMat,}{a matrix representing the normal
  expression}
}
\value{
  the differential expression matrix, it is absolute
  standarized matrix
}
\description{
  DawnNormalize takes a tumor and normal expression matrix
  and returns a standardized absolute differential
  expression matrix. The resulting data is mapped to a standard normal distribution.
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
