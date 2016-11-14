# DawnRank

driver gene identification algorithm DawnRank 2014

## Install:

```r
install.packages("devtools")

library(devtools)

install_github("MartinFXP/DawnRank")

library(DawnRank)
```

## Small scale example:

```r
# load the mutation data
data(brcaExampleMutation)

# load the tumor expression data
data(brcaExampleTumorExpression)

# load the normal expression data
data(brcaExampleNormalExpression)

# load the pathway data
data(brcaExamplePathway)

# load the gold standard
data(goldStandard)

# normalize the tumor and normal data to get the differential expression
normalizedDawn <- DawnNormalize(tumorMat = brcaExampleTumorExpression, normalMat = brcaExampleNormalExpression)

# get the DawnRank Score Get some coffee, this might take a while!
dawnRankScore <- DawnRank(adjMatrix = brcaExamplePathway, mutationMatrix = brcaExampleMutation, 
expressionMatrix = normalizedDawn, mu = 3, goldStandard = goldStandard, parallel = 2)

# look at the DawnRank scores for a few patients
dawnRankFrame <- dawnRankScore[[3]]
head(dawnRankFrame)

# get the aggregate DawnRank scores Get some coffee, this might take a
# while!
aggregateDawnRankScore <- condorcetRanking(scoreMatrix = dawnRankScore[[2]], 
    mutationMatrix = brcaExampleMutation, parallel = 2)

# look at top 10 ranked genes
top10 <- aggregateDawnRankScore[[2]][1:10]
top10

# get the individual cutoff for patient TCGA-A2-A04P
dawnRankFrame$isCGC <- dawnRankFrame$isGoldStandard
library(maxstat)
#NOTE: the latest version of mvnorm should be installed
significance <- patspeccutoff(patient = "TCGA-A2-A04P", ms = dawnRankFrame, 
    default = 95)
# look for signficance. 
significance[[1]]
```

## Reference:

Hou & Ma, 2014

https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0056-8

http://bioen-compbio.bioen.illinois.edu/DawnRank/
