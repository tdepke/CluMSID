% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distancematrix_mds.R
\name{distanceMatrix}
\alias{distanceMatrix}
\title{Create distance matrix from list of spectra}
\usage{
distanceMatrix(speclist, distFun = "cossim", type = c("spectrum",
    "neutral_losses"), mz_tolerance = 1e-05)
}
\arguments{
\item{speclist}{A list of \code{\linkS4class{MS2spectrum}} or
\code{\linkS4class{pseudospectrum}} objects as generated by
\code{\link{extractMS2spectra}} or \code{\link{extractPseudospectra}}.}

\item{distFun}{The distance function to be used. At the moment, only
\code{\link{cossim}} is implemented.}

\item{type}{\code{"spectrum"} (default) for MS2 spectra or MS1 pseudospectra
or \code{"neutral_losses"} for neutral loss patterns.}

\item{mz_tolerance}{The \emph{m/z} tolerance to be used for merging, default
is \code{1e-5}, i.e. +/- 10ppm. If the mass-to-charge ratios of two peaks
differ less than \emph{mz_tolerance}, they are assumed to have the same
\emph{m/z}}
}
\value{
A numeric \code{length(speclist)} by \code{length(speclist)} matrix
    containing pairwise distances (1 - similarity) between all features in
    \code{speclist}. Row and column names are taken from the \code{id} slot
    or, if present, pasted from the \code{id} and \code{annotation} slots of
    the \code{\linkS4class{MS2spectrum}} or
    \code{\linkS4class{pseudospectrum}} objects.
}
\description{
\code{distanceMatrix()} creates a distance matrix from a list of MS2
spectra, MS1 pseudospectra or neutral loss patterns by pairwise comparison
using the specified distance function. This distance matrix is the basis for
CluMSID's data mining functions.
}
\examples{
load(file = system.file("extdata",
    "annotatedSpeclist.RData",
    package = "CluMSIDdata"))

distanceMatrix(annotatedSpeclist[1:20])

}
