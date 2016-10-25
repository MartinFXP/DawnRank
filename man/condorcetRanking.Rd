\name{condorcetRanking}
\alias{condorcetRanking}
\title{condorcetRanking }
\usage{
  condorcetRanking(scoreMatrix, mutationMatrix, pen = 0.85)
}
\arguments{
  \item{scoreMatrix,}{a matrix containing all the given
  DawnRank scores per patient. Rows are genes, columns are
  patients}

  \item{mutationMatrix,}{the mutation matrix. Rows are
  genes, columns are patients}

  \item{pen,}{the penality parameter in the condorcet
  algorithm for missing data. Default 0.85}
}
\value{
  the ranks. A list of 2 including a [[1]] a matrix of all
  pairwise comparisons, [[2]] the final rankings based on
  the Condorcet score
}
\description{
condorcetRanking determines the aggregate rank of a gene over a population.The condorcet rank is used to determine the overall strength of driver gene compared to other drivers by having patients "vote" for a winner using pairwise comparisions.
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
