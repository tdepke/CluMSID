#' Generate a \code{data.frame} with feature information from list of
#' \code{MS2spectrum} objects
#'
#' \code{Featurelist} generates a \code{data.frame} that contains feature ID,
#' precurosur \emph{m/z} and retention time for all features contained in a
#' list of \code{MS2spectrum} objects as produced by \code{extractMS2spectra}
#' and \code{mergeSpecList}. \code{Featurelist} is used internally by
#' \code{\link{writeFeaturelist}}.
#'
#' @param featlist A list of \code{MS2spectrum} objects as produced by
#'   \code{extractMS2spectra} and \code{mergeSpecList}
#'
#' @details Although originally designed for lists of \code{MS2spectrum}
#'   objects, the function also works with lists of \code{pseudospectrum}
#'   objects. In this case, \code{NA} is given for precursor \emph{m/z}.
#'
#' @return A \code{data.frame} that contains feature ID, precurosur \emph{m/z}
#'   (if available) and retention time
#'
#' @examples
#' load(file = system.file("extdata",
#'     "featlist.RData",
#'     package = "CluMSID"))
#'
#' pre_anno <- Featurelist(featlist)
#'
#' @export
Featurelist <- function(featlist){
    id <- c(); mz <- c(); rt <- c()
    for(i in seq_along(featlist)){
        id[i] <- featlist[[i]]@id
        if(methods::.hasSlot(featlist[[i]], "precursor")){
            mz[i] <- featlist[[i]]@precursor
        } else {
            mz[i] <- NA
        }
        rt[i] <- featlist[[i]]@rt
    }
    df <- data.frame(id, mz, rt, stringsAsFactors = FALSE)
    return(df)
}

#' Write feature information from list of \code{MS2spectrum} objects
#'
#' \code{writeFeaturelist} uses \code{\link{Featurelist}} to generate a
#' \code{data.frame} that contains feature ID, precurosur \emph{m/z} and
#' retention time for all features contained in a list of \code{MS2spectrum}
#' objects as produced by \code{extractMS2spectra} and \code{mergeSpecList} and
#' writes it to a csv file.
#'
#' @inheritParams Featurelist
#'
#' @param filename The desired file name of the csv file, default is
#'   \code{"pre_anno.csv"}
#'
#' @details Although originally designed for lists of \code{MS2spectrum}
#'   objects, the function also works with lists of \code{pseudospectrum}
#'   objects. In this case, \code{NA} is given for precursor \emph{m/z}.
#'
#' @return A csv file that contains feature ID, precurosur \emph{m/z} and
#'   retention time. The file has a header but no row names and is separated by
#'   \code{','}.
#'
#' @examples
#' load(file = system.file("extdata",
#'     "featlist.RData",
#'     package = "CluMSID"))
#'
#' writeFeaturelist(featlist, filename = "pre_anno.csv")
#'
#' @export
writeFeaturelist <- function(featlist, filename = "pre_anno.csv"){
    df <- Featurelist(featlist)
    df$annotation <- rep("", length(featlist))
    utils::write.table(
        df,
        file = filename,
        sep = ",",
        row.names = FALSE
    )
}

#' Adding external annotations to list of \code{MS2spectrum} objects
#'
#' \code{addAnnotations} is used to add annotations that have been assigned
#' externally, e.g. by library search, to a list of \code{MS2spectrum} objects
#' as produced by \code{extractMS2spectra} and \code{mergeSpecList}.
#'
#' @inheritParams Featurelist
#'
#' @param annolist A list of annotations, either as a \code{data.frame} or csv
#'   file. The order of features must be the same as in \code{featlist}. Please
#'   see the package vignette for a detailed example!
#'
#' @param annotationColumn The column of \code{annolist} were the annotation is
#'   found. Default is \code{4}, which is the case if
#'   \code{\link{writeFeaturelist}} followed by manual addition of annotations,
#'   e.g. in Excel, is used to generate \code{annolist}.
#'
#' @return A list of \code{MS2spectrum} objects as produced by
#'   \code{extractMS2spectra} and \code{mergeSpecList} with external
#'   annotations added to the \code{annotation} slot of each \code{MS2spectrum}
#'   object.
#'
#' @examples
#' load(file = system.file("extdata",
#'     "featlist.RData",
#'     package = "CluMSID"))
#'
#' addAnnotations(featlist, system.file("extdata",
#'                 "post_anno.csv",
#'                 package = "CluMSID"),
#'                 annotationColumn = 4)
#'
#' @export
addAnnotations <- function(featlist, annolist, annotationColumn = 4){
    if(is.data.frame(annolist)){
        ident <- annolist
    } else {
        ident <-
            utils::read.csv(file = annolist, stringsAsFactors = FALSE)
    }
    stopifnot(length(featlist) == nrow(annolist))
    for(i in seq_along(featlist)){
        if(ident[i, annotationColumn] != ""){
            featlist[[i]]@annotation <- ident[i, annotationColumn]
        }
    }
    return(featlist)
}
