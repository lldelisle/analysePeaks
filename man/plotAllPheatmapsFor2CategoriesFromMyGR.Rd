% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingFunctions.R
\name{plotAllPheatmapsFor2CategoriesFromMyGR}
\alias{plotAllPheatmapsFor2CategoriesFromMyGR}
\title{Plot pheatmaps for categories of 2 specific features for a Set and a Ref}
\usage{
plotAllPheatmapsFor2CategoriesFromMyGR(
  myGRs,
  nameOfColWithCate1,
  cateNames1,
  what1,
  nameOfColWithCate2,
  cateNames2,
  what2,
  stringSet,
  nameOfRef,
  fontsize = 10,
  display_numbers = T,
  whichPheatmaps = 1:5,
  plotProportion = F
)
}
\arguments{
\item{myGRs}{a GRangeList of 2 elements which are the set and the ref and should have 2 meta columns `nameOfColWithCate1` and `nameOfColWithCate2`}

\item{nameOfColWithCate1}{a string which gives the column to display}

\item{cateNames1}{the names of the categories, in the good order, for the first feature}

\item{what1}{the name of the first measure}

\item{nameOfColWithCate2}{a string which gives the column to display for the first feature}

\item{cateNames2}{the names of the categories, in the good order, for the second feature}

\item{what2}{the name of the second measure}

\item{stringSet}{the name of the set for the second feature}

\item{nameOfRef}{the name of the reference}

\item{fontsize}{base fontsize for the heatmaps (default is 10)}

\item{display_numbers}{logical whether to display the numbers in the pheatmaps (default is TRUE)}

\item{whichPheatmaps}{a vector containing the pheatmaps to plot between 1 to 5 with:
1=shared by all Set and at least one Ref or shared by all Set and Ref
2=shared by all Set none of the Ref or shared by all Set no Ref
3=shared by all Set: all or all shared by all Set
4=shared by all Ref in not a single of the set or Ref in not a single of the set
5=shared by all Ref: all or all Ref
(default is 1:5)}

\item{plotProportion}{logical whether to display the pheatmap 1, 2, and 4 should be plot as proportion of pheatmap 3 and 5 (default is FALSE)}
}
\value{
Plot but do not return anything
}
\description{
Plot pheatmaps for categories of 2 specific features for a Set and a Ref
}
