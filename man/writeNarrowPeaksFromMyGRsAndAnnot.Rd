% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otherFunctions.R
\name{writeNarrowPeaksFromMyGRsAndAnnot}
\alias{writeNarrowPeaksFromMyGRsAndAnnot}
\title{Extract common and specific peaks in narrowPeak format}
\usage{
writeNarrowPeaksFromMyGRsAndAnnot(
  myGRAndAttributes,
  outputDirectory,
  useOriginalNarrowPeak = F
)
}
\arguments{
\item{myGRAndAttributes}{Should be the output of \link[analysePeaks]{createMyGRs} or \link[analysePeaks]{createMyGRsUsingSimpleOverlap}}

\item{outputDirectory}{a path to a folder where will be written the narrowPeaks files.}

\item{useOriginalNarrowPeak}{logical to specify if instead of the coordinates in myGRAndAttributes which are by default the summit +/-150bp. You want the original coordinates as in the input (default is F).}
}
\value{
Does not return anything but write 4 files into the `outputDirectory`:
Ref_in_at_least_one_Set.narrowPeak or shared_by_all_Ref_and_at_least_one_Set.narrowPeak
Ref_in_not_a_single_Set.narrowPeak or share_by_all_Ref_and_not_in_a_single_Set.narrowPeak
share_by_all_Set_no_Ref.narrowPeak or share_by_all_Set_and_not_in_a_single_Ref.narrowPeak
shared_by_all_Set_and_Ref.narrowPeak or shared_by_all_Set_and_at_least_one_Ref.narrowPeak
}
\description{
Extract common and specific peaks in narrowPeak format
}
