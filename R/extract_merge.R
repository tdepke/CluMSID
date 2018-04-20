#' Merge spectra with m/z tolerance
#'
#' \code{mergeTolerance()} merges two spectra by identifying common peaks with a
#' given m/z tolerance. It can be used with \code{Reduce()} to merge more than
#' two spectra.
#'
#' @param x,y MS2 spectra as objects of class \code{matrix} with m/z in the
#'   first column and intensity in the second.
#'
#' @param tolerance The m/z tolerance used for merging. If two peaks are within
#'   tolerance, they are regarded as the same. Defaults to \code{1e-5}, i.e.
#'   10ppm.
#'
#' @return A matrix with m/z in the first column and separate columns for
#'   intensities in the respective spectra. If peaks were merged, their m/z
#'   corresponds to the mean of the two original m/z.
#'
#' @keywords internal
mergeTolerance <- function(x, y, tolerance = 1e-5) {
  mrg <- merge(x, y, by = "V1", all = T)
  mrg[is.na(mrg)] <- 0
  i <- 1
  while (!is.na(mrg[(i + 1), 1])) {
    if (abs(mrg[i, 1] - mrg[(i + 1), 1]) <= mrg[i, 1] * tolerance) {
      mrg[i, 1] <- (mrg[i, 1] + mrg[(i + 1), 1]) / 2
      mrg[i,-1] <- mrg[i,-1] + mrg[(i + 1),-1]
      mrg <- mrg[-(i + 1),]
      i <- i + 1
      colnames(mrg) <-
        c("V1", 2:ncol(mrg)) #suppresses error warning 'duplicate column names'
    } else {
      i <- i + 1
    }
  }
  mrg
}

#' Extract MS2 spectra from raw data files
#'
#' \code{extractMS2spectra()} is used to extract MS2 spectra from raw data
#' files, e.g. mzXML files.
#'
#' @param MSfile An LC-MS/MS raw data file in one of the non-proprietary formats
#'   that can be parsed by \code{mzR}, e.g. mzXML or mzML.
#'
#' @param min_peaks Minimum number of peaks in MS2 spectrum, defaults to
#'   \code{2}. Spectra with less than \code{min_peaks} fragment peaks will be
#'   ignored and not extracted.
#'
#' @param recalibrate_precursor Logical, defaults to \code{FALSE}. Applicable
#'   only for files that were exported to mzXML using a deprecated version of
#'   Bruker Compass Xport. If set to \code{TRUE}, the precursor m/z will be
#'   recalculated from the respective fragment m/z in the MS2 spectrum. For
#'   details, see Depke et al. 2017.
#'
#' @param RTlims Retention time interval for the extraction of spectra. Provide
#'   as numeric vector of length 2. Spectra with retention time <
#'   \code{RTlims[1]} or > \code{RTlims[2]} will be ignored.
#'
#' @return A \code{list} with objects of class \code{MS2spectrum}, containing
#'   MS2 spectra extracted from the raw data.
#'
#' @export
extractMS2spectra <- function(MSfile, min_peaks = 2,
                              recalibrate_precursor = FALSE,
                              RTlims = NULL){
  aa <- mzR::openMSfile(MSfile, backend = "Ramp")

  mslvl <- c()
  for (z in 1:length(aa)) {
    mslvl[z] <- mzR::header(aa, z)$msLevel
  }
  try(if(length(mslvl[mslvl == 2]) < 1) stop("The file does not contain MS2 spectra."))

  spectra <- list()
  for (z in 1:length(aa)) {
    spectra[[z]] <- mzR::peaks(aa, z)
  }
  ms2log <- mslvl == 2
  ms2spectra <- spectra[ms2log]

  vec <- c()
  for (k in 1:length(ms2spectra)) {
    vec[k] <- (nrow(ms2spectra[[k]]) >= min_peaks)
  }
  ms2spectra2 <- ms2spectra[vec]

  pmz <- mzR::header(aa)$precursorMZ

  if(recalibrate_precursor == TRUE){
    new.pmz <- 0
    for (i in 2:length(pmz)) {
      if (pmz[i] == 0) {
        x <-
          0
      } else {
        if (pmz[(i - 1)] == 0) {
          x <-
            mzR::peaks(aa, (i - 1))[which.min(abs(pmz[i] - mzR::peaks(aa, (i - 1))[, 1])), 1]
        } else {
          if (pmz[(i - 2)] == 0) {
            x <-
              mzR::peaks(aa, (i - 2))[which.min(abs(pmz[i] - mzR::peaks(aa, (i - 2))[, 1])), 1]
          } else {
            if (pmz[(i - 3)] == 0) {
              x <-
                mzR::peaks(aa, (i - 3))[which.min(abs(pmz[i] - mzR::peaks(aa, (i - 3))[, 1])), 1]
            } else {
              x <- NA
            }
          }
        }
      }
      if (x == 0 ||
          ((abs(x - pmz[i]) / pmz[i]) * 1e06) <= 100) {
        new.pmz[i] <- x
      } else {
        new.pmz[i] <- NA
      }
    }
  } else {
    new.pmz <- pmz
  }


  precursor <- cbind(new.pmz, mzR::header(aa)$retentionTime)
  precursor2 <- precursor[ms2log,][vec,]

  if(!is.null(RTlims)){
    cutrt <- precursor2[, 2] < (RTlims[2] * 60) & precursor2[, 2] > (RTlims[1] * 60)
    precursormzrt <- precursor2[cutrt,]
    ms2list <- ms2spectra2[cutrt]
  } else {
    precursormzrt <- precursor2
    ms2list <- ms2spectra2
  }

  output <- list()
  for(e in 1:length(ms2list)){
    output[[e]] <- methods::new("MS2spectrum",
                       precursor = as.numeric(precursormzrt[e,1]),
                       rt = as.numeric(precursormzrt[e,2]),
                       spectrum = ms2list[[e]])
  }
  return(output)
  close(aa)
}

#' Merge list of spectra
#'
#' \code{mergeSpecList()} is an accessory function used only inside
#' \code{mergeMS2spectra}.
#'
#' @param speclist A \code{list} of \code{MS2spectrum} objects to be merged.
#'
#' @param tolerance The m/z tolerance to be used for merging.
#'
#' @return A \code{list} of the same length as \code{speclist} containing merged
#'   spectra as \code{MS2spectrum} objects. If multiple spectra contribute to
#'   one consensus spectrum, than this consensus spectrum is contained in the
#'   list multiple times at the respective positions of the contributing
#'   spectra.
#'
#' @keywords internal
mergeSpecList <- function(speclist, tolerance = 1e-5) {
  mergeToleranceX <- function(x,y){
    mergeTolerance(x,y,tolerance = tolerance)
  } #to circumvent the problem that Reduce() cannot handle additional arguments
  mrgls <- list()
  ident <- c()
  for (s in 1:length(speclist)) {
    ident[s] <- speclist[[s]]@id
  }
  for (z in 1:length(speclist)) {
    z0 <- c()
    for (j in 1:(z - 1)) {
      z0[j] <- ident[z] == ident[j]
    }
    if (z != 1 & any(z0)) {
      mrgls[[z]] <- mrgls[[which(z0 == T)[1]]]
    } else {
      if (sum(ident == ident[z]) > 1) {
        zl <- speclist[ident == ident[z]]
        for(d in 1:length(zl)){
          zl[[d]] <- zl[[d]]@spectrum
        }
        z1 <- as.matrix(Reduce(mergeToleranceX, zl)) #Problem: cannot include "tolerance" arg in Reduce
        z1[is.na(z1)] <- 0
        z2 <-
          cbind(z1[, 1], round((rowSums(z1) - z1[, 1]) / ncol(z1[,-1])))
        dimnames(z2) <- NULL
        mrgls[[z]] <- methods::new("MS2spectrum",
                          id = ident[z],
                          precursor = as.numeric(speclist[[z]]@precursor),
                          rt = as.numeric(speclist[[z]]@rt),
                          spectrum = z2)
      } else {
        mrgls[[z]] <- speclist[[z]]
      }
    }
  }
  mrgls
}

#' Generate neutral loss patterns from MS2 spectra
#'
#' \code{neutrallossPatterns} generates neutral loss patterns from MS2 spectra
#' and adds them to \code{\linkS4class{MS2spectrum}} objects in the slot
#' \code{neutral_losses}.
#'
#' @param x an object of class \code{\linkS4class{MS2spectrum}} that contains an
#'   MS2 spectrum in the \code{spectrum} slot
#'
#' @return an object of class \code{\linkS4class{MS2spectrum}} with a neutral
#'   loss pattern in the \code{neutral_losses} slot
#'
#' @keywords internal
neutrallossPatterns <- function(x){
  temp.nl <- cbind((x@precursor - x@spectrum[, 1]),
                   x@spectrum[, 2])
  temp.nl <- subset(temp.nl, temp.nl[, 1] >= (x@precursor * 1e-5))
  #include unfragmented precursor??
  x@neutral_losses <- temp.nl
  return(x)
}

#' Merge MS2 spectra with or without external peak table
#'
#' \code{mergeMS2spectra} is used to merge MS2 spectra that come from the same precursor.
#' It does so either by grouping spectra of the same precursor \emph{m/z} that fall into
#' a defined retention time window (\code{rt_tolerance}) or by grouping spectra to peaks
#' from an externally supplied peak table. Please see the vignette for more details.
#'
#' @param ms2list A \code{list} of \code{MS2spectrum} objects to be merged.
#'
#' @param mz_tolerance The \emph{m/z} tolerance to be used for merging, default is \code{1e-5},
#' i.e. +/- 10ppm. If the mass-to-charge ratios of two peaks differ less than
#' \emph{mz_tolerance}, they are assumed to have the same \emph{m/z}
#'
#' @param rt_tolerance The retention time tolerance used for merging features. If used without
#' a peak table, \code{rt_tolerance} is the maximum retention time difference between to subsequent
#' spectra of the same precursor \emph{m/z} with which they are still assumed to belong to the same feature
#' If used with an external peak table, \code{rt_tolerance} is the maximum retention time difference between
#' a spectrum and a peak in the peak table with which the spectrum is still considered to belong to that peak.
#'
#' @param peaktable An external peak table, e.g. from XCMS, that serves as a template for grouping spectra.
#' The peaktable must be a three-column \code{data.frame} with feature ID, \emph{m/z} and retention time for each
#' peak/feature.
#'
#' @param exclude_unmatched If an external peak table is used: Should spectra that do not match
#' to any peak/feature in the peak table be exclude from the resulting list?
#'
#' @return A merged list of \code{\linkS4class{MS2spectrum}} objects.
#'
#' @export
mergeMS2spectra <- function(ms2list,
                            mz_tolerance = 1e-5,
                            rt_tolerance = 30,
                            peaktable = NULL,
                            exclude_unmatched = FALSE){
  flist <- list()
  mz <- c()
  for(k in 1:length(ms2list)){
    mz[k] <- ms2list[[k]]@precursor
  }
  rt <- c()
  for(k in 1:length(ms2list)){
    rt[k] <- ms2list[[k]]@rt
  }

  mz1 <- cbind(mz, rt)
  if(any(is.na(mz1))){
    stop("NAs in either mz or rt slot in at least one object!")
  }

  #if no sample table is provided, the original
  #algorithm is used to summarise spectra/features
  if(is.null(peaktable)){
    while (nrow(mz1) >= 1) {
      l1 <- abs(mz1[1, 1] - mz1[, 1]) <= mz1[1, 1] * mz_tolerance
      l2 <- matrix(mz1[c(l1, l1)], ncol = 2)
      l3 <- diff(l2[, 2])
      l4 <- c(0, which(l3 > 30), nrow(l2))
      l5 <- list()
      for (i in 1:(length(l4) - 1)) {
        l5[[i]] <- l2[(l4[i] + 1):(l4[i + 1]),]
      }
      flist <- append(flist, l5)
      mz1 <- matrix(mz1[c(!l1, !l1)], ncol = 2)
    }


    for (i in 1:length(flist)) {
      if (is.matrix(flist[[i]])) {
        flist[[i]] <- cbind(flist[[i]],
                            rep(stats::median(flist[[i]][, 1]), times = nrow(flist[[i]])),
                            rep(stats::median(flist[[i]][, 2]), times = nrow(flist[[i]])))
      } else {
        flist[[i]] <- c(flist[[i]], flist[[i]])
      }
    }
    medmzrt <- c()
    for (i in 1:length(flist)) {
      medmzrt <- rbind(medmzrt, flist[[i]])
    }
    medmzrt <-
      as.data.frame(medmzrt)
    colnames(medmzrt) <- c("mz", "rt", "med.mz", "med.rt")

    medmzrt$id <- paste("M", round((medmzrt$med.mz), 2),
                        "T", round((medmzrt$med.rt), 2), sep = "")

    medmzrt <- medmzrt[order(medmzrt$rt),]

    for (g in 1:length(ms2list)){
      temp <- which((ms2list[[g]]@precursor == medmzrt$mz)
                    & (ms2list[[g]]@rt == medmzrt$rt))
      ms2list[[g]]@id <- medmzrt$id[temp]
      ms2list[[g]]@precursor <- round(medmzrt$med.mz[temp], 4)
      ms2list[[g]]@rt <- round(medmzrt$med.rt[temp], 2)
    }


  } else { #when a peak table is used ...


    matr <- matrix(data = NA, ncol = 3, nrow = nrow(mz1))
    for(e in 1:nrow(mz1)){
      l1 <- abs(mz1[e, 1] - peaktable[,2]) <= mz1[e, 1] * mz_tolerance &
        abs(mz1[e, 2] - peaktable[, 3]) <= rt_tolerance
      if(sum(l1) == 0){
        matr[e,] <- c(paste0("no_match_", e), mz1[e, 1], mz1[e, 2])
      } else if(sum(l1) == 1){
        matr[e,] <- unlist(peaktable[l1,])
      } else if(sum(l1) > 1){
        matr[e,] <- unlist(peaktable[l1,][which.min(unlist(abs(peaktable[l1,3] - mz1[e, 2]))),])
      }
    }

    mz2 <- matr[grepl(pattern = "no_match_", x = matr[,1]),2:3]
    mz2 <- matrix(data = sapply(mz2, as.numeric), ncol = 2)

    while (nrow(mz2) >= 1) {
      l1 <- abs(mz2[1, 1] - mz2[, 1]) <= mz2[1, 1] * mz_tolerance
      l2 <- matrix(mz2[c(l1, l1)], ncol = 2)
      l3 <- diff(l2[, 2])
      l4 <- c(0, which(l3 > 30), nrow(l2))
      l5 <- list()
      for (i in 1:(length(l4) - 1)) {
        l5[[i]] <- l2[(l4[i] + 1):(l4[i + 1]),]
      }
      flist <- append(flist, l5)
      mz2 <- matrix(mz2[c(!l1, !l1)], ncol = 2)
    }

    for (i in 1:length(flist)) {
      if (is.matrix(flist[[i]])) {
        flist[[i]] <- cbind(flist[[i]],
                            rep(stats::median(flist[[i]][, 1]), times = nrow(flist[[i]])),
                            rep(stats::median(flist[[i]][, 2]), times = nrow(flist[[i]])))
      } else {
        flist[[i]] <- c(flist[[i]], flist[[i]])
      }
    }
    medmzrt <- c()
    for (i in 1:length(flist)) {
      medmzrt <- rbind(medmzrt, flist[[i]])
    }
    medmzrt <-
      as.data.frame(medmzrt)
    colnames(medmzrt) <- c("mz", "rt", "med.mz", "med.rt")

    medmzrt$id <- paste("xM", round((medmzrt$med.mz), 2),
                        "T", round((medmzrt$med.rt), 2), sep = "")

    matr2 <- as.matrix(medmzrt)[,c(5,3,4)]

    matr[grepl(pattern = "no_match_", x = matr[,1]),] <- matr2

    for (g in 1:length(ms2list)){
      ms2list[[g]]@id <- matr[g,1]
      ms2list[[g]]@precursor <- round(as.numeric(matr[g,2]), 4)
      ms2list[[g]]@rt <- round(as.numeric(matr[g,3]), 2)
    }
  }

  mergedlist <- mergeSpecList(ms2list, tolerance = mz_tolerance)
  shortlist <- mergedlist[!duplicated(mergedlist)]

  for(u in 1:length(shortlist)){
    shortlist[[u]] <- neutrallossPatterns(shortlist[[u]])
  }

  if(exclude_unmatched == TRUE){
    w <- 1
    while(w <= length(shortlist)){
      if(grepl(pattern = "xM", x = shortlist[[w]]@id)){
        shortlist[w] <- NULL
      } else w <- w+1
    }
  }

  return(shortlist)
}
