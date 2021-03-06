% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_merge.R
\name{extractMS2spectra}
\alias{extractMS2spectra}
\title{Extract MS2 spectra from raw data files}
\usage{
extractMS2spectra(MSfile, min_peaks = 2, recalibrate_precursor = FALSE,
    RTlims = NULL)
}
\arguments{
\item{MSfile}{An LC-MS/MS raw data file in one of the non-proprietary
formats that can be parsed by \code{mzR}, e.g. mzXML or mzML.}

\item{min_peaks}{Minimum number of peaks in MS2 spectrum, defaults to
\code{2}. Spectra with less than \code{min_peaks} fragment peaks will be
ignored and not extracted.}

\item{recalibrate_precursor}{Logical, defaults to \code{FALSE}. Applicable
only for files that were exported to mzXML using a deprecated version of
Bruker Compass Xport (< 3.0.13). If set to \code{TRUE}, the precursor m/z
will be recalculated from the respective fragment m/z in the MS2 spectrum.
For details, see Depke et al. 2017.}

\item{RTlims}{Retention time interval for the extraction of spectra. Provide
as numeric vector of length 2. Spectra with retention time <
\code{RTlims[1]} or > \code{RTlims[2]} will be ignored.}
}
\value{
A \code{list} with objects of class \code{MS2spectrum}, containing
    MS2 spectra extracted from the raw data.
}
\description{
\code{extractMS2spectra()} is used to extract MS2 spectra from raw data
files, e.g. mzXML files.
}
\examples{
my_spectra <- extractMS2spectra(MSfile = system.file("extdata",
                                "PoolA_R_SE.mzXML",
                                package = "CluMSIDdata"),
                                min_peaks = 4, RTlims = c(0,10))

}
