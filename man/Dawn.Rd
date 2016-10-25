\name{Dawn}
\alias{Dawn}
\title{Dawn}
\usage{
  Dawn(adjMatrix, expressionVector, mutationVector,
    damping, maxit = 100, epsilon = 1e-04,
    goldStandard = NULL, patientTag = "defaultPatient")
}
\arguments{
  \item{adjMatrix,}{the adjacency}

  \item{expressionVector,}{the normalized expression
  vector}

  \item{mutationVector,}{a logical vector containing
  mutation information}

  \item{damping,}{the damping vector}

  \item{maxit,}{the maximum number of iterations to use,
  default 100}

  \item{epsilon,}{the lower magnitude cutoff, default
  0.0001}

  \item{goldStandard,}{A list of common driver genes, used
  as a comparison. This is optional, default=NULL}

  \item{patientTag,}{an index for the patients}
}
\value{
  the ranks. A list of 3 including a [[1]] output of all
  the ranks, [[2]] mutated ranks, [[3]] the steps of
  convergence
}
\description{
 The shell of the DawnRank algorithm for single patient. Dawn determines rankings of mutated genes for only one patient. Dawn is called in a loop by DawnRank.
 }

