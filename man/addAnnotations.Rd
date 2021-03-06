% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{addAnnotations}
\alias{addAnnotations}
\title{Adding external annotations to list of \code{MS2spectrum} objects}
\usage{
addAnnotations(featlist, annolist, annotationColumn = 4)
}
\arguments{
\item{featlist}{A list of \code{MS2spectrum} objects as produced by
\code{extractMS2spectra} and \code{mergeSpecList}}

\item{annolist}{A list of annotations, either as a \code{data.frame} or csv
file. The order of features must be the same as in \code{featlist}. Please
see the package vignette for a detailed example!}

\item{annotationColumn}{The column of \code{annolist} were the annotation is
found. Default is \code{4}, which is the case if
\code{\link{writeFeaturelist}} followed by manual addition of annotations,
e.g. in Excel, is used to generate \code{annolist}.}
}
\value{
A list of \code{MS2spectrum} objects as produced by
    \code{extractMS2spectra} and \code{mergeSpecList} with external
    annotations added to the \code{annotation} slot of each \code{MS2spectrum}
    object.
}
\description{
\code{addAnnotations} is used to add annotations that have been assigned
externally, e.g. by library search, to a list of \code{MS2spectrum} objects
as produced by \code{extractMS2spectra} and \code{mergeSpecList}.
}
\examples{
load(file = system.file("extdata",
    "featlist.RData",
    package = "CluMSIDdata"))

addAnnotations(featlist, system.file("extdata",
                "post_anno.csv",
                package = "CluMSIDdata"),
                annotationColumn = 4)

}
