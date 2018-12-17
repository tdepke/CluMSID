#' Accessor functions for individual slots of
#' \code{\linkS4class{MS2spectrum}} and
#' \code{\linkS4class{pseudospectrum}} objects
#'
#' @param x An object of class \code{\linkS4class{MS2spectrum}}
#'   or \code{\linkS4class{pseudospectrum}}
#'
#' @return The value of the respective slot of the object
#'   (\code{id}, \code{annotation}, \code{precursor},
#'   \code{rt}, \code{spectrum}, \code{neutral_losses})
#'
#' @name accessors
NULL

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' accessID(annotatedSpeclist[[1]])
#'
#' @export
accessID <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@id)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' accessAnnotation(annotatedSpeclist[[1]])
#'
#' @export
accessAnnotation <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@annotation)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' accessPrecursor(annotatedSpeclist[[1]])
#'
#' @export
accessPrecursor <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@precursor)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'      "annotatedSpeclist.RData",
#'      package = "CluMSIDdata"))
#'
#' accessRT(annotatedSpeclist[[1]])
#'
#' @export
accessRT <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@rt)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'      "annotatedSpeclist.RData",
#'      package = "CluMSIDdata"))
#'
#' accessPolarity(annotatedSpeclist[[1]])
#'
#' @export
accessPolarity <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@polarity)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' accessSpectrum(annotatedSpeclist[[1]])
#'
#' @export
accessSpectrum <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@spectrum)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' accessNeutralLosses(annotatedSpeclist[[1]])
#'
#' @export
accessNeutralLosses <- function(x){
    stopifnot(is(x, "MS2spectrum") | is(x, "pseudospectrum"))
    return(x@neutral_losses)
}
