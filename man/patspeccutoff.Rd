\name{patspeccutoff}
\alias{patspeccutoff}
\title{patspeccuttoff}
\usage{
  patspeccutoff(patient, ms, default = 0.95)
}
\arguments{
  \item{patient,}{the patient of interest}

  \item{ms,}{a dataframe containing dawn rank output
  (specifically, the third element of DawnRank output)}

  \item{default,}{the default cutoff in case there are not
  enough common drivers to make a clear conclusion}
}
\value{
  a data frame containing information regarding the genetic
  information of a specific patient, and the cutoff itself
}
\description{
  patspeccuttoff determines the score cutoffs to distinguish drivers and passenger mutations of a single patient using Maximially-Selected Rank Statistics. These cutoffs determine a threshold for each patient that maximizes the number of known genes that score above that threshold. 
}
\references{
Torsten Hothorn and Berthold Lausen. On maximally selected rank statistics. R News, 2(1):3–5, 2002.
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
