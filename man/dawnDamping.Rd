\name{dawnDamping}
\alias{dawnDamping}
\title{dawnDamping}
\usage{
  dawnDamping(adjMatrix, mu)
}
\arguments{
  \item{adjMatrix,}{the adjacency matrix}

  \item{mu,}{the free parameter to fit}
}
\value{
  the vector for individualized damping factors
}
\description{
  dawnDamping calculates the damping factor used in DawnRank using the DirichletRank dynamic damping factor.DawnDamping can be called by DawnRank should the user supply the mu parameter.
}

