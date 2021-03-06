% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{networkplot}
\alias{networkplot}
\title{Correlation network from spectral similarity data}
\usage{
networkplot(distmat, interactive = FALSE, show_labels = FALSE,
    label_size = 1.5, highlight_annotated = FALSE,
    min_similarity = 0.1, exclude_singletons = FALSE)
}
\arguments{
\item{distmat}{A distance matrix as generated by
\code{\link{distanceMatrix}}.}

\item{interactive}{Logical, defaults to \code{FALSE}. If \code{TRUE}, an
interactive visualisation is generated using \pkg{plotly}.}

\item{show_labels}{Logical, defaults to \code{FALSE}. If \code{TRUE},
feature IDs are printed as labels in the network plot. Argument has no
effect if \code{interactive} is \code{TRUE} (because in this case, labels
are displayed on mouse-over).}

\item{label_size}{Numeric, defaults to \code{1.5}. If \code{show_labels} is
\code{TRUE} and \code{interactive} is \code{FALSE}, \code{label_size}
defines the size of labels in the plot.}

\item{highlight_annotated}{Logical, defaults to \code{FALSE}. If
\code{TRUE}, points for features for which an annotation was added before
using \code{\link{distanceMatrix}} are highlighted by red colour, while
other points are grey in the network plot.}

\item{min_similarity}{Numeric, defaults to \code{0.1}. The minimum spectral
contrast angle (see\code{\link{cossim}}) that is considered a spectral
similarity and hence a connection in the network.}

\item{exclude_singletons}{Logical, defaults to \code{FALSE}. If \code{TRUE},
features that have no connection to any other feature will not be
displayed in the network plot.}
}
\value{
A network plot generated with the help of
    \code{\link[network]{network}}, \code{\link[GGally]{ggnet2}} and, if
    interactive, \code{\link[plotly]{ggplotly}}. Edge weights correspond to
    spectral similarities.
}
\description{
\code{networkplot()} is used to generate correlation networks from
spectral similarity data. An interactive visualisation can be produced using
\pkg{plotly}.
}
\examples{
load(file = system.file("extdata",
    "distmat.RData",
    package = "CluMSIDdata"))

networkplot(distmat[1:50,1:50], show_labels = TRUE,
                exclude_singletons = TRUE)

}
