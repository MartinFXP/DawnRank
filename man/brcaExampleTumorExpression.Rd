\name{brcaExampleTumorExpression}
\alias{brcaExampleTumorExpression}
\docType{data}
\title{
brcaExampleTumorExpression
}
\description{
Tumor Expression data for an example Breast Cancer Dataset from the TCGA.The data is presented in matrix with genes being the rows and patients being columns. This dataset contains the unnormalized tumor expression.
}
\usage{data(brcaExampleTumorExpression)}
\format{
  The format is:
 num [1:1492, 1:533] 0.0578 0.3625 0.7008 0.3775 0.4605 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
  ..$ : chr [1:533] "TCGA-A1-A0SD" "TCGA-A1-A0SE" "TCGA-A1-A0SH" "TCGA-A1-A0SJ" ...
}
\details{
Expression data for an example Breast Cancer Dataset. Rows are genes, columns are patients
}
\source{
https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
}
\references{
The Cancer Genome Atlas. Comprehensive molecular portraits of human breast tumors. 
Nature. 490 (7418):61-70.}
\examples{
data(brcaExampleTumorExpression)
## maybe str(brcaExampleTumorExpression) ;
plot(brcaExampleTumorExpression)
}
\keyword{datasets}
