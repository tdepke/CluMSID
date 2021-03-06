#' Match one spectrum against a set of spectra
#'
#' \code{getSimilarities} calculates the similarities of one spectrum or
#' neutral loss pattern to a set of other spectra or neutral loss patterns.
#'
#' @param spec The spectrum to be compared to other spectra. Can be either an
#'   object of class \code{\linkS4class{MS2spectrum}} or a two-column numerical
#'   matrix that contains fragment mass-to-charge ratios in the first and
#'   intensities in the second column.
#'
#' @param speclist The set of spectra to which \code{spec} is to be compared.
#'   Must be a list where every entry is an object of class
#'   \code{\linkS4class{MS2spectrum}}. Can be generated from an mzXML file with
#'   \code{\link{extractMS2spectra}} and \code{\link{mergeMS2spectra}} or
#'   constructed using \code{new("MS2spectrum", ...)} for every list entry (see
#'   vignette for details).
#'
#' @param type Specifies whether MS2 spectra or neutral loss patterns are to be
#'   compared. Must be either \code{'spectrum'} (default) or
#'   \code{'neutral_losses'}.
#'
#' @param hits_only Logical that indicates whether the result should contain
#'   only similarities greater than zero.
#'
#' @return A named vector with similarities of \code{spec} to all spectra or
#'   neutral loss patterns in \code{speclist}.
#'
#' @importFrom methods new is
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#' getSimilarities(annotatedSpeclist[[137]],
#'                 annotatedSpeclist, hits_only = TRUE)
#'
#' @export
getSimilarities <- function(spec,
                            speclist,
                            type = c("spectrum", "neutral_losses"),
                            hits_only = FALSE){
    stopifnot((is(spec, "MS2spectrum") |
                is(spec, "pseudospectrum") |
                is.matrix(spec)) & is.list(speclist) &
                all(vapply(X = speclist,
                FUN = function(x) is(x, "MS2spectrum")|is(x, "pseudospectrum"),
                FUN.VALUE = logical(1L))))
    type <- match.arg(type)
    if(is.matrix(spec)){
        if(type == "spectrum"){
            spec <- methods::new("MS2spectrum", spectrum = spec)
        } else {
            spec <- methods::new("MS2spectrum", neutral_losses = spec)
        }
    }
    simvec <- c()
    for(k in seq_along(speclist)){
        simvec[k] <- cossim(spec, speclist[[k]], type = type)
        names(simvec)[k] <- speclist[[k]]@id
    }
    if(hits_only){
        simvec <- simvec[simvec > 0 & !is.na(simvec)]
    }
    return(simvec)
}

#' Find spectra that contain a specific fragment
#'
#' \code{findFragment} is used to find spectra that contain a specific fragment
#' ion. Its sister function is \code{\link{findNL}}, which finds specific
#' neutral losses. Both functions work analogous to \code{\link{getSpectrum}}.
#'
#' @param featlist a list that contains only objects of class
#'   \code{\linkS4class{MS2spectrum}}
#'
#' @param mz The mass-to-charge ratio of the fragment ion of interest.
#'
#' @param tolerance The \emph{m/z} tolerance for the fragment ion search.
#'   Default is \code{1E-05}, i.e. +/- 10ppm.
#'
#' @return If the respective fragment is only found in one spectrum, the output
#'   is an object of class \code{\linkS4class{MS2spectrum}}; if it is found in
#'   more than one spectrum, the output is a list of
#'   \code{\linkS4class{MS2spectrum}} objects.
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#' putativeAQs <- findFragment(annotatedSpeclist, 159.068)
#'
#' @export
findFragment <- function(featlist, mz, tolerance = 1E-05){
    subsetter <- vapply(X = featlist,
                        FUN = function(x) any(abs(accessSpectrum(x)[,1]
                                                - mz) <= mz * tolerance),
                        FUN.VALUE = logical(1L))
    message(paste(sum(subsetter),
                    "spectra were found that contain a fragment of m/z",
                    mz, "+/-", tolerance * 1E06, "ppm."))
    return(featlist[subsetter])
}

#' Find spectra that contain a specific neutral loss
#'
#' \code{findNL} is used to find spectra that contain a specific neutral loss.
#' Its sister function is \code{\link{findFragment}}, which finds specific
#' fragment ions. Both functions work analogous to \code{\link{getSpectrum}}.
#'
#' @param featlist a list that contains only objects of class
#'   \code{\linkS4class{MS2spectrum}}
#'
#' @param mz The mass-to-charge ratio of the neutral loss of interest.
#'
#' @param tolerance The \emph{m/z} tolerance for the neutral loss search.
#'   Default is \code{1E-05}, i.e. +/- 10ppm.
#'
#' @return If the respective neutral loss is only found in one spectrum, the
#'   output is an object of class \code{\linkS4class{MS2spectrum}}; if it is
#'   found in more than one spectrum, the output is a list of
#'   \code{\linkS4class{MS2spectrum}} objects.
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#' findNL(annotatedSpeclist, 212.009)
#'
#' @export
findNL <- function(featlist, mz, tolerance = 1E-05){
    subsetter <- vapply(X = featlist,
                        FUN = function(x) any(abs(accessNeutralLosses(x)[,1]
                                                    - mz) <= mz * tolerance),
                        FUN.VALUE = logical(1L))
    message(paste(sum(subsetter),
                    "neutral loss patterns were found that",
                    "contain a neutral loss of m/z",
                    mz, "+/-", tolerance * 1E06, "ppm."))
    return(featlist[subsetter])
}

#' Access individual spectra from a list of spectra by various slot entries
#'
#' As accessing S4 objects within lists is not trivial, \code{getSpectrum} can
#' be used to access individual or several \code{\linkS4class{MS2spectrum}}
#' objects by their slot entries.
#'
#' @param featlist a list that contains only objects of class
#'   \code{\linkS4class{MS2spectrum}}
#'
#' @param slot The slot to be searched (invalid slot arguments will produce
#'   errors). Possible values are: \itemize{ \item \code{'id'} \item
#'   \code{'annotation'} \item \code{'precursor'} (\emph{m/z} of precursor ion)
#'   \item \code{'rt'} (retention time of precursor) }
#'
#' @param what the search term or number, must be character for \code{'id'} and
#'   \code{'annotation'} and numeric for \code{'precursor'} and \code{'rt'} See
#'   vignette for examples.
#'
#' @param mz.tol the tolerance used for precursor ion *m/z* searches, defaults
#'   to \code{1E-05} (+/- 10ppm)
#'
#' @param rt.tol the tolerance used for precursor ion retention time searches,
#'   defaults to 30s; high values can be used to specify retention time ranges
#'   (see vignette for example)
#'
#' @return If the only one spectrum matches the search criteria, the output is
#'   an object of class \code{\linkS4class{MS2spectrum}}; if more than one
#'   spectrum matches, the output is a list of \code{\linkS4class{MS2spectrum}}
#'   objects.
#'
#' @importFrom methods slot
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' getSpectrum(annotatedSpeclist, "annotation", "pyocyanin")
#'
#' getSpectrum(annotatedSpeclist, "id", "M244.17T796.4")
#'
#' getSpectrum(annotatedSpeclist, "precursor", 286.18, mz.tol = 1E-03)
#'
#' six_eight <- getSpectrum(annotatedSpeclist, "rt", 420, rt.tol = 60)
#'
#' @export
getSpectrum <- function(featlist, slot, what, mz.tol = 1E-05, rt.tol = 30){
    stopifnot(slot %in% c("id", "annotation", "precursor", "rt"))
    qexp <- switch(
        slot,
        "id" = quote(what %in% methods::slot(e, slot)),
        "annotation" = quote(what %in% methods::slot(e, slot)),
        "precursor" = quote(abs(what - methods::slot(e, slot)) <= mz.tol),
        "rt" = quote(abs(what - methods::slot(e, slot)) <= rt.tol)
    )
    subsetter <- vapply(X = featlist,
                        FUN = function(e) eval(qexp),
                        FUN.VALUE = logical(1L))
    if(sum(subsetter) == 0){
        message("No spectrum with that ", slot, ".", sep = "")
    } else if(sum(subsetter) == 1){
        return(featlist[subsetter][[1]]) # return as MS2spectrum object
    } else {
        return(featlist[subsetter]) # return as list
    }
}

#' Separate spectra with different polarities from the same run
#'
#' Using \code{splitPolarities}, spectra with different polarities from the same
#' run can be separated, e.g. when processing spectra recorded with
#' polarity-switching.
#'
#' @param ms2list A list of \code{\linkS4class{MS2spectrum}} objects as produced
#'   by \code{\link{extractMS2spectra}}.
#'
#' @param polarity The polarity of spectra to be analysed, must be
#'   \code{"positive"} or \code{"negative"}.
#'
#' @return A list of \code{\linkS4class{MS2spectrum}} objects that contains only
#'   spectra with the given \code{polarity}.
#'
#' @examples
#' my_spectra <- extractMS2spectra(MSfile = system.file("extdata",
#'                                 "PoolA_R_SE.mzXML",
#'                                 package = "CluMSIDdata"),
#'                                 min_peaks = 4, RTlims = c(0,5))
#'
#' my_positive_spectra <- splitPolarities(my_spectra, "positive")
#'
#' @export
splitPolarities <- function(ms2list, polarity = c("positive", "negative")){
    stopifnot(polarity %in% c("positive", "negative"))
    subvec <- vapply(FUN = function (x) accessPolarity(x) == polarity,
                        X = ms2list,
                        FUN.VALUE = logical(1L))
    return(ms2list[subvec])
}

#' Create a basic plot of MS2 spectra
#'
#' \code{specplot} creates a very basic plot of MS2 spectra from
#' \code{\linkS4class{MS2spectrum}} or \code{\linkS4class{pseudospectrum}}
#' objects.
#'
#' @param spec An object of class \code{\linkS4class{MS2spectrum}} or
#'   \code{\linkS4class{pseudospectrum}}
#'
#' @param ... Additional graphical parameters to be passed to \code{plot()}
#'
#' @return A plot of the MS2 spectrum saved in the \code{spectrum} slot of
#'   \code{spec}.
#'
#' @importFrom graphics plot text
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' specplot(annotatedSpeclist[[1]])
#'
#' @export
specplot <- function(spec, ...) {
    stopifnot(is(spec, "MS2spectrum") | is(spec, "pseudospectrum"))
    params <- list(...)
    if(!("main" %in% names(params))) {
        params$main <- paste("id:", spec@id, " - ", "rt:", spec@rt)
    }
    if(!("sub" %in% names(params))) params$sub <- spec@annotation
    params$x <- spec@spectrum[, 1]
    params$y <- spec@spectrum[, 2] / max(spec@spectrum[, 2])
    params$type <- "h"
    params$xlim <- c(0, (max(spec@spectrum[, 1]) * 1.1))
    params$xaxs <- "i"
    params$xlab <- expression(italic(m / z))
    params$ylim <- c(0, 1.1)
    params$yaxs <- "i"
    params$ylab <- "intensity relative to base peak"
    do.call(graphics::plot, params)
    params$x <- (spec@spectrum[, 1])[(spec@spectrum[, 2] /
                                        max(spec@spectrum[, 2])) > 0.1]
    params$y <- (spec@spectrum[, 2] /
                    max(spec@spectrum[, 2]))[
                        (spec@spectrum[, 2] /
                            max(spec@spectrum[, 2])) > 0.1]
    if(!("labels" %in% names(params))) {
        params$labels <- round((spec@spectrum[, 1])
                                [(spec@spectrum[, 2] /
                                    max(spec@spectrum[, 2])) > 0.1], 4)
    }
    if(!("pos" %in% names(params))) params$pos <- 3
    if(!("cex" %in% names(params))) params$cex <- 0.75
    params$type <- NULL
    do.call(graphics::text, params)
}
