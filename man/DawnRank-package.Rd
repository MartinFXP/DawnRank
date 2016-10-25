\name{DawnRank-package}
\alias{DawnRank-package}
\alias{DawnRank}
\docType{package}
\title{
DawnRank
}
\description{
DawnRank calculates the importance of potential driver mutations in individual patientes using a pathway approach inspired by the PageRank algorithm. Mutations in genes that are more highly connected and connect to differentially expressed downstream genes are more likely to be driver mutations. DawnRank can calculate the importance of candidate driver mutations in a single patient as well as make predictions on the common drivers from a patient population using Condorcet Rank Aggregation.
}
\details{
\tabular{ll}{
Package: \tab DawnRank\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-10-15\cr
License: \tab GPL\cr
}
}
\author{
Jack Hou <jackhou2@illinois.edu>, Jian Ma

}
\references{
}
\keyword{ package }
\seealso{
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
