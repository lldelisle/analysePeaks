% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otherFunctions.R
\name{createMyGRsUsingSimpleOverlap}
\alias{createMyGRsUsingSimpleOverlap}
\title{Create a GRList with a reference GR and a sample GR annotated with GRanges provided in arguments. Here this is not a one-to-one match between the Set and the Reference but just an overlap.}
\usage{
createMyGRsUsingSimpleOverlap(
  mySet,
  allExpe,
  ovCluster = NULL,
  grLScores = NULL,
  grLDistance = NULL,
  myRef = c("ChIP"),
  useSummitPMFlanking = T,
  flankingSize = 150
)
}
\arguments{
\item{mySet}{a vector with the name of the samples to compare to the ref}

\item{allExpe}{a GRangeList with at least each sample in `mySet` and `myRef`}

\item{ovCluster}{a dataframe summurizing the cluster by overlap if previously calculated (default is NULL)}

\item{grLScores}{a GRangeList with annotations for which you would like to know the best score overlapping the GRanges of experiments}

\item{grLDistance}{a GRangeList with annotations for which you would like to know the distance to the closest item}

\item{myRef}{among the samples which one is the reference, the others will be merged as "replicates"}

\item{useSummitPMFlanking}{logical value to specify if you prefer to use the full region of the peak or only the summit +/- the `flankingSize`}

\item{flankingSize}{integer value used to extend the summit if `useSummitPMFlanking` is TRUE}
}
\value{
A list with:
`stringSet` a string with all samples other than the reference separated by comma
`myGRs` a GRangeList of 2 elements which are
      1. The peaks which are shared by all samples of mySet
      2. The peaks which are shared by all samples of myRef
`nameOfRef` the name of the reference used
`namesOfGrLScores` the names of grLScores
`namesOfGrLDistance` the names of grLDistance
`useSummitPMFlanking` the value of useSummitPMFlanking which was used
`flankingSize` the value of flankingSize which was used
}
\description{
Create a GRList with a reference GR and a sample GR annotated with GRanges provided in arguments. Here this is not a one-to-one match between the Set and the Reference but just an overlap.
}
\details{
The GRanges of myGRs will have in mcols a new information which is
bestNameOfTheGRScore for each GR in grLScores and
distanceToTheNearestNameOfTheGR for each GR in grLDistance
The first item of myGRs have the coordinates of the first item of `mySet` and will also have inNoneOfTheRef.
The second item of myGRs have the coordinates of the first item of `myRef` and will also have inNoneOfTheSet.
Every overlap was considered and used to make "cluster", only the peak with the higest score whithin the same cluster is used in the output.
}
