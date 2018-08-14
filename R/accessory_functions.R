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
#' @export
getSimilarities <- function(spec,
                            speclist,
                            type = "spectrum",
                            hits_only = FALSE){
  if(!(type %in% c("spectrum", "neutral_losses"))){
    stop("'type' must be either 'spectrum' (default) or 'neutral_losses'!")
  }
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
  if(hits_only == TRUE){
    simvec <- simvec[simvec > 0]
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
#' @export
findFragment <- function(featlist, mz, tolerance = 1E-05){
  subsetter <- c()
  for(i in seq_along(featlist)){
    m <- featlist[[i]]@spectrum
    subsetter[i] <- any(abs(m[,1] - mz) <= mz * tolerance)
  }
  message(cat(sum(subsetter),
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
#' @export
findNL <- function(featlist, mz, tolerance = 1E-05){
  subsetter <- c()
  for(i in seq_along(featlist)){
    m <- featlist[[i]]@neutral_losses
    subsetter[i] <- any(abs(m[,1] - mz) <= mz * tolerance)
  }
  message(cat(sum(subsetter),
              "neutral loss patterns were found that
              contain a neutral loss of m/z",
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
#' @export
getSpectrum <- function(featlist, slot, what, mz.tol = 1E-05, rt.tol = 30){
  subsetter <- c()

  if(slot %in% c("id", "annotation")){
    for(i in seq_along(featlist)){
      m <- methods::slot(featlist[[i]], slot)
      subsetter[i] <- what %in% m
    }
  } else if(slot == "precursor"){
    for(i in seq_along(featlist)){
      m <- methods::slot(featlist[[i]], slot)
      subsetter[i] <- abs(what - m) <= mz.tol
    }
  } else if(slot == "rt"){
    for(i in seq_along(featlist)){
      m <- methods::slot(featlist[[i]], slot)
      subsetter[i] <- abs(what - m) <= rt.tol
    }
  } else stop("invalid slot selected")


  if(sum(subsetter) == 0){
    cat("No spectrum with that ", slot, ".", sep = "")
  } else if(sum(subsetter) == 1){
    return(featlist[subsetter][[1]])
  } else {
    return(featlist[subsetter])
  }
}
