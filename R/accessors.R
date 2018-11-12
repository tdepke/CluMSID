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
#' access_id(annotatedSpeclist[[1]])
#'
#' @export
access_id <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@id)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' access_annotation(annotatedSpeclist[[1]])
#'
#' @export
access_annotation <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@annotation)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' access_precursor(annotatedSpeclist[[1]])
#'
#' @export
access_precursor <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@precursor)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'      "annotatedSpeclist.RData",
#'      package = "CluMSIDdata"))
#'
#' access_rt(annotatedSpeclist[[1]])
#'
#' @export
access_rt <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@rt)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'      "annotatedSpeclist.RData",
#'      package = "CluMSIDdata"))
#'
#' access_polarity(annotatedSpeclist[[1]])
#'
#' @export
access_polarity <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@polarity)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' access_spectrum(annotatedSpeclist[[1]])
#'
#' @export
access_spectrum <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@spectrum)
}

#' @rdname accessors
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' access_neutral_losses(annotatedSpeclist[[1]])
#'
#' @export
access_neutral_losses <- function(x){
    stopifnot(class(x) %in% c("MS2spectrum", "pseudospectrum"))
    return(x@neutral_losses)
}
