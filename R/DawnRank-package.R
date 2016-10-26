

#' brcaExampleMutation
#' 
#' Mutation data for an example Breast Cancer Dataset
#' 
#' Mutation data for an example Breast Cancer Dataset. Rows are genes, columns
#' are patients
#' 
#' @name brcaExampleMutation
#' @docType data
#' @format The format is: num [1:1492, 1:533] 0 0 0 0 0 0 0 0 0 0 ...  -
#' attr(*, "dimnames")=List of 2 ..$ : chr [1:1492] "PML" "FGF3" "HLA-A"
#' "CACNA1F" ...  ..$ : chr [1:533] "TCGA-A1-A0SD" "TCGA-A1-A0SE"
#' "TCGA-A1-A0SH" "TCGA-A1-A0SJ" ...
#' @references TCGA
#' @source https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
#' @keywords datasets
#' @examples
#' 
#' data(brcaExampleMutation)
#' ## maybe str(brcaExampleMutation) ; plot(brcaExampleMutation) ...
#' 
NULL





#' brcaExampleNormalExpression
#' 
#' Normal Expression data for an example Breast Cancer Dataset
#' 
#' Expression data for an example Breast Cancer Dataset. Rows are genes,
#' columns are patients
#' 
#' @name brcaExampleNormalExpression
#' @docType data
#' @format The format is: num [1:1492, 1:62] 0.24 0.122 0.56 -0.231 1.189 ...
#' - attr(*, "dimnames")=List of 2 ..$ : chr [1:1492] "PML" "FGF3" "HLA-A"
#' "CACNA1F" ...  ..$ : chr [1:62] "TCGA-A7-A0CE" "TCGA-A7-A0CH" "TCGA-A7-A0D9"
#' "TCGA-A7-A0DB" ...
#' @references TCGA
#' @source https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
#' @keywords datasets
#' @examples
#' 
#' data(brcaExampleNormalExpression)
#' ## maybe str(brcaExampleNormalExpression) ; plot(brcaExampleNormalExpression) ...
#' 
NULL





#' brcaExamplePathway
#' 
#' Pathway data for an example Breast Cancer Dataset
#' 
#' An adjacency matrix detailing interactions of an example pathway
#' 
#' @name brcaExamplePathway
#' @docType data
#' @format The format is: num [1:1492, 1:1492] 0 0 0 0 0 0 0 0 0 0 ...  -
#' attr(*, "dimnames")=List of 2 ..$ : chr [1:1492] "PML" "FGF3" "HLA-A"
#' "CACNA1F" ...  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
#' @references KEGG
#' @source http://www.genome.jp/kegg/
#' @keywords datasets
#' @examples
#' 
#' data(brcaExamplePathway)
#' ## maybe str(brcaExamplePathway) ; plot(brcaExamplePathway) ...
#' 
NULL





#' brcaExampleTumorExpression
#' 
#' Tumor Expression data for an example Breast Cancer Dataset
#' 
#' Expression data for an example Breast Cancer Dataset. Rows are genes,
#' columns are patients
#' 
#' @name brcaExampleTumorExpression
#' @docType data
#' @format The format is: num [1:1492, 1:533] 0.0578 0.3625 0.7008 0.3775
#' 0.4605 ...  - attr(*, "dimnames")=List of 2 ..$ : chr [1:1492] "PML" "FGF3"
#' "HLA-A" "CACNA1F" ...  ..$ : chr [1:533] "TCGA-A1-A0SD" "TCGA-A1-A0SE"
#' "TCGA-A1-A0SH" "TCGA-A1-A0SJ" ...
#' @references TCGA
#' @source https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
#' @keywords datasets
#' @examples
#' 
#' data(brcaExampleTumorExpression)
#' ## maybe str(brcaExampleTumorExpression) ; plot(brcaExampleTumorExpression) ...
#' 
NULL





#' DawnRank
#' 
#' DawnRank calculates the pathway impact a mutated gene has on a individual
#' patient using a PageRank inspired algorithm. Mutations in genes with greater
#' pathway impact are likely to be driver genes in cancer. DawnRank can
#' calculate individual pathway impacts and predict drivers from groups of
#' patients
#' 
#' \tabular{ll}{ Package: \tab DawnRank\cr Type: \tab Package\cr Version: \tab
#' 1.1\cr Date: \tab 2013-10-15\cr License: \tab GPL\cr }
#' 
#' @name DawnRank-package
#' @aliases DawnRank-package DawnRank
#' @docType package
#' @author Jack Hou <jackhou2@@illinois.edu> Maintainer: Martin Pirkl <martin.pirkl@bsse.ethz.ch>
#' @seealso
#' @references
#' @keywords package
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
NULL





#' goldStandard
#' 
#' The Cancer Gene Census Gold standard for common driver genes.
#' 
#' The Cancer Gene Census Gold standard for common driver genes.This is a
#' vector
#' 
#' @name goldStandard
#' @docType data
#' @format The format is: chr [1:147] "ABL1" "ABL2" "ACSL3" "AKAP9" "AKT1"
#' "AKT2" "APC" "ARID1A" ...
#' @references COSMIC
#' @source http://cancer.sanger.ac.uk/cancergenome/projects/census/
#' @keywords datasets
#' @examples
#' 
#' data(goldStandard)
#' ## maybe str(goldStandard) ; plot(goldStandard) ...
#' 
NULL



