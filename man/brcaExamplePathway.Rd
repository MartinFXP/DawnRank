\name{brcaExamplePathway}
\alias{brcaExamplePathway}
\docType{data}
\title{
brcaExamplePathway
}
\description{
Pathway data (from KEGG) for an example Breast Cancer Dataset representing all networks in the to the KEGG "Pathways in Cancer". The pathway is in the form of an adjacency matrix describing curated gene-gene interactions.
}
\usage{data(brcaExamplePathway)}
\format{
  The format is:
 num [1:1492, 1:1492] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
  ..$ : chr [1:1492] "PML" "FGF3" "HLA-A" "CACNA1F" ...
}
\details{
An adjacency matrix detailing interactions of an example pathway
}
\source{
http://www.genome.jp/kegg/
}
\references{
Kanehisa, M., Goto, S., Sato, Y., Furumichi, M., and Tanabe, M.; KEGG for integration and interpretation of large-scale molecular datasets. Nucleic Acids Res. 40, D109-D114 (2012).
}
\examples{
data(brcaExamplePathway)
## maybe str(brcaExamplePathway) ; plot(brcaExamplePathway)...
}
\keyword{datasets}
