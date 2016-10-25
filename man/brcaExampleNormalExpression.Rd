\name{brcaExampleNormalExpression}
\alias{brcaExampleNormalExpression}
\docType{data}
\title{
brcaExampleNormalExpression
}
\description{
Normal Expression data for an example Breast Cancer Dataset from TCGA. The data is in matrix with genes being the rows and patients being columns. This dataset contains the unnormalized normal expression of normal BRCA TCGA samples.
}
\usage{data(brcaExampleNormalExpression)}
\format{
  The format is:
 num [1:1492, 1:62] 0.24 0.122 0.56 -0.231 1.189 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
  ..$ : chr [1:62] "TCGA-A7-A0CE" "TCGA-A7-A0CH" "TCGA-A7-A0D9" "TCGA-A7-A0DB" ...
}
\details{
Expression data for an example Breast Cancer Dataset. Rows are genes, columns are patients
}
\source{
https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
}
\references{
The Cancer Genome Atlas. Comprehensive molecular portraits of human breast tumors. 
Nature. 490 (7418):61-70.
}
\examples{
data(brcaExampleNormalExpression)
## maybe str(brcaExampleNormalExpression) ; 
}
\keyword{datasets}
