#' Extract pseudospectra
#'
#' \code{extractPseudospectra()} is used to extract MS1 pseudospectra from
#' \pkg{CAMERA} output.
#'
#' @param x \pkg{CAMERA} output that contains information on pseudospectra.
#' Can either be of class \code{data.frame} or
#' \code{\link[CAMERA]{xsAnnotate}}.
#' It is recommended to use either \code{xsAnnotate} objects or
#' \code{data.frame}s generated from XCMSonline results tables but
#' other \code{data.frame}s are possible.
#'
#' @param min_peaks Minimum number of peaks in pseudospectrum, defaults to
#' \code{1}. See \code{\link{extractMS2spectra}}.
#'
#' @param intensity_columns Numeric, defaults to \code{NULL}.
#' If a \code{data.frame} is used as input which has not been
#' generated from an XCMSonline results table, the indices
#' of the columns that contain the peak intensities in the
#' different samples have to be indicated as
#' \code{intensity_columns}.
#'
#' @return A list of pseudospectra, stored as objects of class
#' \code{\linkS4class{pseudospectrum}}, analogous to the output of
#' \code{\link{extractMS2spectra}}.
#'
#' @examples
#' pstable <- readr::read_delim(file = system.file("extdata",
#'                                 "TD035_XCMS.annotated.diffreport.tsv",
#'                                 package = "CluMSID"), delim = "\t")
#'
#' pseudospeclist <- extractPseudospectra(pstable, min_peaks = 2)
#'
#' @export
extractPseudospectra <- function(x, min_peaks = 1, intensity_columns = NULL){
    ##different actions depending on class of x
    if(is.data.frame(x)){

        ##extract variables
        pcg <- x$pcgroup

        #mz and rt can have different names ...
        if("mz" %in% colnames(x)){
            mz <- x$mz
        } else if("mzmed" %in% colnames(x)){
            mz <- x$mzmed
        } else stop("The column name for m/z must be either 'mz' or 'mzmed'!")

        if("rt" %in% colnames(x)){
            rt <- x$rt
        } else if("rtmed" %in% colnames(x)){
            rt <- x$rtmed
        } else stop("The column name for retention time
                must be either 'rt' or 'rtmed'!")

        ##XCMSonline output has maximum intensity as a column,
        ##this is easiest to use for intensity
        ##otherwise, intensity columns have to be indicated!
        if("maxint" %in% colnames(x)){
            maxint <- x$maxint
        } else if(!is.null(intensity_columns)){
            maxint <- c()
            for(i in seq_len(nrow(x))){
                maxint[i] <- max(x[i, intensity_columns], na.rm = TRUE)
            }
        } else stop("If x does not contain a 'maxint' column,
                please indicate the indices of the columns
                that contain intensity values!")

        ##Create actual list of pseudospectra
        pseudospeclist <- list()

        for(i in unique(pcg)){
            pseudospeclist[[i]] <- methods::new("pseudospectrum",
                                                id = i,
                                                rt = mean(rt[pcg == i]),
                                                spectrum =
                                                    cbind(mz[pcg == i],
                                                            maxint[pcg == i]))
        }

        ##If filtered data are used, some pcgroups can be empty,
        ##those are discarded in this step
        names(pseudospeclist) <- seq_along(pseudospeclist)
        pseudospeclist <- Filter(Negate(is.null), pseudospeclist)

    } else if(class(x) == "xsAnnotate"){

        ##xsAnnotate objects have a different structure depending on whether
        ##they were generated from single or multiple mzXML files
        if("maxo" %in% colnames(x@groupInfo)){

            pseudospeclist <- list()
            for(i in seq_along(x@pspectra)){

                spc <- cbind(x@groupInfo[x@pspectra[[i]],"mz"],
                                x@groupInfo[x@pspectra[[i]],"maxo"])

                pseudospeclist[[i]] <- methods::new("pseudospectrum",
                                                    id = i,
                                                    rt =
                                                        stats::median(
                                                        x@groupInfo[
                                                            x@pspectra[[i]],
                                                            "rt"]),
                                                    spectrum = spc)
            }
        } else {

            ##automatically find intensity columns by
            ##fuzzy matching of sample names
            temp <- list()
            for(j in seq_along(rownames(x@xcmsSet@phenoData))){
                temp[[j]] <- agrep(rownames(x@xcmsSet@phenoData)[j],
                                    colnames(x@groupInfo))
            }
            sv <- unique(unlist(temp))

            ##as rowMax() doesn't have a na.rm feature,
            ##we need to get rid of all NAs
            x@groupInfo[,sv][is.na(x@groupInfo[,sv])] <- 0

            pseudospeclist <- list()
            for(i in seq_along(x@pspectra)){
                if(length(x@pspectra[[i]])>1){
                    spc <- cbind(x@groupInfo[x@pspectra[[i]],"mz"],
                                    Biobase::rowMax(
                                        x@groupInfo[x@pspectra[[i]],sv]))
                } else {
                    spc <- cbind(x@groupInfo[x@pspectra[[i]],"mz"],
                                    max(x@groupInfo[x@pspectra[[i]],sv],
                                        na.rm = TRUE))
                }
                pseudospeclist[[i]] <- methods::new("pseudospectrum",
                                                    id = i,
                                                    rt =
                                                    stats::median(
                                                    x@groupInfo[x@
                                                                pspectra[[i]],
                                                                "rt"]),
                                                    spectrum = spc)
            }

        }

    } else stop("'x' must be either a data.frame
                or an object of class 'xsAnnotate'!")

    ##Possibility to filter pseudospectra by minimum number of peaks
    psvec <- c()
    for(i in seq_along(pseudospeclist)) {
        psvec[i] <- nrow(pseudospeclist[[i]]@spectrum) > min_peaks
    }

    return(pseudospeclist[psvec])
}
