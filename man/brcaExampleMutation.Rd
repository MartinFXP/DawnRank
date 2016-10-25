\name{brcaExampleMutation}
\alias{brcaExampleMutation}
\docType{data}
\title{
brcaExampleMutation
}
\description{
Mutation data for a Breast Cancer Dataset from The Cancer Gene Atlas (TCGA). The data is presented in matrix with genes being the rows and patients being columns. A binary 0 or 1 is used to indicate whether or not a gene is mutated in that patient.
}
\usage{data(brcaExampleMutation)}
\format{
  The format is:
 num [1:1492, 1:533] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
  ..$ : chr [1:533] "TCGA-A1-A0SD" "TCGA-A1-A0SE" "TCGA-A1-A0SH" "TCGA-A1-A0SJ" ...
}
\details{
Mutation data for an example Breast Cancer Dataset. Rows are genes, columns are patients
}
\source{
https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
}
\references{
The Cancer Genome Atlas. Comprehensive molecular portraits of human breast tumors. 
Nature. 490 (7418):61-70.
}
\examples{
data(brcaExampleMutation)
## maybe str(brcaExampleMutation) ; plot(brcaExampleMutation) ...
}
\keyword{datasets}
