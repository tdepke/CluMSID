#' A custom S4 class for MS2 spectra, neutral loss patterns and
#' respective metainformation
#'
#' @slot id a character string similar to the ID used by XCMSonline
#'   or the ID given in a predefined peak list
#'
#' @slot annotation a character string containing a user-defined
#'   annotation, defaults to empty
#'
#' @slot precursor (median) \emph{m/z} of the spectrum's precursor ion
#'
#' @slot rt (median) retention time of the spectrum's precursor ion
#'
#' @slot polarity the ionisation polarity,
#'    \code{"positive"} or \code{"negative"}
#'
#' @slot spectrum the actual MS2 spectrum as two-column matrix
#'   (column 1 is (median) \emph{m/z}, column 2 is (median) intensity of the
#'   product ions)
#'
#' @slot neutral_losses a neutral loss pattern generated by subtracting the
#'   product ion mass-to-charge ratios from the precursor \emph{m/z} in a
#'   matrix format analogous to the \code{spectrum} slot
#'
#' @importFrom methods setClass
#'
#' @export
setClass("MS2spectrum",
        slots = list(   id = "character",
                        annotation = "character",
                        precursor = "numeric",
                        rt = "numeric",
                        polarity = "character",
                        spectrum = "matrix",
                        neutral_losses = "matrix"))

#' A custom S4 class for MS1 pseudospectra and respective
#' metainformation
#'
#' @slot id a the \code{"pcgroup"} number assigned by \pkg{CAMERA}
#'
#' @slot annotation a character string containing a user-defined
#'   annotation, defaults to empty
#'
#' @slot rt (median) retention time of the ions contained in the
#'   pseudospectrum
#'
#' @slot spectrum the actual MS1 pseudospectrum as two-column
#'   matrix (column 1 is (median) \emph{m/z}, column 2 is (median)
#'   intensity of the ions)
#'
#' @importFrom methods setClass
#'
#' @export
setClass("pseudospectrum",
        slots = list(   id = "numeric",
                        annotation = "character",
                        rt = "numeric",
                        spectrum = "matrix"))

#' @describeIn MS2spectrum A show generic for \code{MS2spectra}.
#'
#' @param object An object of class \code{\linkS4class{MS2spectrum}}
#'
#' @return Prints information from the object slots with exception of
#'   'spectrum' and 'neutral_losses' where only a summary is given.
#'
#' @importFrom methods setMethod show
#'
#' @exportMethod show
setMethod("show",
            "MS2spectrum",
            function(object){
                cat('An object of class "MS2spectrum"', '\n',
                    'id:', object@id, '\n',
                    'annotation:', object@annotation, '\n',
                    'precursor:', format(object@precursor, nsmall = 4), '\n',
                    'retention time:', object@rt, '\n',
                    'polarity:', object@polarity, '\n',
                    'MS2 spectrum with', nrow(object@spectrum),
                    'fragment peaks',
                    '\n', 'neutral loss pattern with',
                    nrow(object@neutral_losses),
                    'neutral losses'
                )
            })

#' Convert spectra from \pkg{MSnbase} classes
#'
#' @param x An object of class \code{\link[MSnbase:Spectrum-class]{Spectrum}}
#'   or \code{\link[MSnbase:Spectrum-class]{Spectrum2}}
#'
#' @return An object of class \code{\linkS4class{MS2spectrum}}
#'
#' @importFrom methods new
#' @importClassesFrom MSnbase Spectrum Spectrum2
#'
#' @examples
#' #Load a "Spectrum2" object from MSnbase
#' library(MSnbase)
#' sp <- itraqdata[["X1"]]
#' #Convert this object to "MS2spectrum" class
#' new_sp <- as.MS2spectrum(sp)
#' #Or alternatively:
#' new_sp <- as(sp, "MS2spectrum")
#'
#' @export
as.MS2spectrum <- function(x) as(x, "MS2spectrum")
setAs("Spectrum", "MS2spectrum",
        function(from) {
            methods::new(
                "MS2spectrum",
                precursor = from@precursorMz,
                rt = from@rt,
                spectrum = cbind(from@mz, from@intensity)
            )
        })
setAs("Spectrum2", "MS2spectrum",
        function(from) {
            methods::new(
                "MS2spectrum",
                precursor = from@precursorMz,
                rt = from@rt,
                spectrum = cbind(from@mz, from@intensity)
            )
        })

#' Calculate cosine similarity between two spectra
#'
#' \code{cossim()} calculates the cosine of the spectral constrast angle as a
#' measure for the similarity of two spectra.
#'
#' @param x,y MS2 spectra, either as \code{matrix},
#'   \code{\linkS4class{MS2spectrum}} or \code{\linkS4class{pseudospectrum}}
#'   objects. \code{x} and \code{y} must have the same class.
#'
#' @param type Whether similarity between spectra (\code{"spectrum"}, default)
#'   or neutral loss patterns (\code{"neutral_losses"}) is to be compared
#'
#' @param mzTolerance The m/z tolerance used for merging. If two fragment peaks
#'   are within tolerance, they are regarded as the same. Defaults to
#'   \code{1e-5}, i.e. 10ppm.
#'
#' @return The cosine similarity of \code{x} and \code{y}
#'
#' @examples
#' load(file = system.file("extdata",
#'     "annotatedSpeclist.RData",
#'     package = "CluMSIDdata"))
#'
#' cossim(annotatedSpeclist[[1]], annotatedSpeclist[[2]])
#'
#' @export
cossim <- function(x, y, type = c("spectrum", "neutral_losses"),
                    mzTolerance = 1e-5) {
    colnames(x) <- NULL
    colnames(y) <- NULL
    mm <- mergeTolerance(x, y, tolerance = mzTolerance)
    sum(sqrt(mm[, 2]) * sqrt(mm[, 3])) /
        (sqrt(sum(mm[, 2])) * sqrt(sum(mm[, 3])))
}
setGeneric("cossim")

#' @describeIn cossim \code{cossim} method for
#'   \code{\linkS4class{MS2spectrum}} objects
#'
#' @importFrom methods setMethod
#'
#' @exportMethod cossim
setMethod("cossim",
            c(x = "MS2spectrum", y = "MS2spectrum"),
            function(x, y, type, mzTolerance){
                type <- match.arg(type)
                if(type == "spectrum"){
                    cossim(x@spectrum, y@spectrum, type = type,
                            mzTolerance = mzTolerance)
                } else if(type == "neutral_losses"){
                    cossim(x@neutral_losses, y@neutral_losses,
                            type = type, mzTolerance = mzTolerance)
                } else stop("'type' must be either 'spectrum' (default)
                            or 'neutral_losses'")
            })

#' @describeIn cossim \code{cossim} method for
#'   \code{\linkS4class{pseudospectrum}} objects
#'
#' @importFrom methods setMethod
#'
#' @exportMethod cossim
setMethod(  "cossim",
            c(x = "pseudospectrum", y = "pseudospectrum"),
            function(x, y, type, mzTolerance){
                type <- match.arg(type)
                if(type == "spectrum"){
                    cossim(x@spectrum, y@spectrum,
                            type = type, mzTolerance = mzTolerance)
                } else stop("with pseudospectra, 'type' must be 'spectrum'!")
            })

#' @describeIn MS2spectrum Method for\code{MSnbase::precursorMz}
#'   for \code{\linkS4class{MS2spectrum}} objects. Accesses \code{precursor}
#'   slot and returns precursor \emph{m/z} as a numeric.
#'
#' @importFrom MSnbase precursorMz
#'
#' @exportMethod precursorMz
setMethod("precursorMz",
            "MS2spectrum",
            function(object) {
                return(object@precursor)
            })

#' @describeIn MS2spectrum Method for\code{MSnbase::rtime}
#'   for \code{\linkS4class{MS2spectrum}} objects. Accesses \code{rt}
#'   slot and returns retention time as a numeric.
#'
#' @importFrom MSnbase rtime
#'
#' @exportMethod rtime
setMethod("rtime",
            "MS2spectrum",
            function(object) {
                return(object@rt)
            })

#' @describeIn MS2spectrum Method for\code{MSnbase::intensity}
#'   for \code{\linkS4class{MS2spectrum}} objects. Accesses \code{spectrum}
#'   slot and returns the intensity column as a numeric vector.
#'
#' @importFrom MSnbase intensity
#'
#' @exportMethod intensity
setMethod("intensity",
            "MS2spectrum",
            function(object) {
                return(object@spectrum[,2])
            })

#' @describeIn MS2spectrum Method for\code{MSnbase::mz}
#'   for \code{\linkS4class{MS2spectrum}} objects. Accesses \code{spectrum}
#'   slot and returns the \emph{m/z} column as a numeric vector.
#'
#' @importFrom MSnbase mz
#'
#' @exportMethod mz
setMethod("mz",
            "MS2spectrum",
            function(object) {
                return(object@spectrum[,1])
            })

#' @describeIn MS2spectrum Method for\code{MSnbase::mz}
#'   for \code{\linkS4class{MS2spectrum}} objects. Accesses \code{spectrum}
#'   slot and returns the number of peaks as a numeric.
#'
#' @importFrom MSnbase peaksCount
#'
#' @exportMethod peaksCount
setMethod("peaksCount",
            "MS2spectrum",
            function(object) {
                return(nrow(object@spectrum))
            })
