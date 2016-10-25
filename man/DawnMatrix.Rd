\name{DawnMatrix}
\alias{DawnMatrix}
\title{DawnMatrix}
\usage{
  DawnMatrix(adjMatrix)
}
\arguments{
  \item{adjMatrix,}{the adjacency matrix}
}
\value{
  the transition matrix
}
\description{
  DawnMatrix creates a transition matrix for DawnRank where the transition probability of each node is determined by the number of incoming edges to that node. 
  
  }
\references{
Page, Lawrence and Brin, Sergey and Motwani, Rajeev and Winograd, Terry (1999) The PageRank Citation Ranking: Bringing Order to the Web. Technical Report.
}
