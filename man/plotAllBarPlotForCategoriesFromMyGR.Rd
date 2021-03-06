% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingFunctions.R
\name{plotAllBarPlotForCategoriesFromMyGR}
\alias{plotAllBarPlotForCategoriesFromMyGR}
\title{Plot barplot for categories of a specific feature for a Set and a Ref}
\usage{
plotAllBarPlotForCategoriesFromMyGR(
  myGRs,
  nameOfColWithCate,
  cateNames,
  stringSet,
  what,
  nameOfRef,
  only5categoriesBP = FALSE
)
}
\arguments{
\item{myGRs}{a GRangeList of 2 elements which are the set and the ref and should have as meta column `nameOfColWithCate`}

\item{nameOfColWithCate}{a string which gives the column to display}

\item{cateNames}{the names of the categories, in the good order.}

\item{stringSet}{the name of the set}

\item{what}{the name of the measure}

\item{nameOfRef}{the name of the reference}

\item{only5categoriesBP}{logical whether to plot only the barplots with 5 categories (default is FALSE).}
}
\value{
Plot but do not return anything
}
\description{
Plot barplot for categories of a specific feature for a Set and a Ref
}
