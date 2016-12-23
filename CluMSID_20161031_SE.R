# CluMSID 0.1

library(mzR)
library(iontree)
library(ape)
library(dbscan)
library(RColorBrewer)

## Open connection to mzXML file
aa <- openMSfile("PoolA_R_SE.mzXML")

## Count MS2 spectra (QC)
mslvl <- c()
for (z in 1:length(aa)) {
  mslvl[z] <- header(aa, z)$msLevel
}
length(mslvl[mslvl == 2])

## Extract MS2 spectra
spectra <- list()
for (z in 1:length(aa)) {
  spectra[[z]] <- peaks(aa, z)
}
ms2log <- mslvl == 2
ms2spectra <- spectra[ms2log]

## Create list with spectra containing 2 or more peaks
vec <- c()
for (k in 1:length(ms2spectra)) {
  vec[k] <- (nrow(ms2spectra[[k]]) >= 2)
}
ms2spectra2 <- ms2spectra[vec]

## Correct uncalibrated precursor masses
pmz <- header(aa)$precursorMZ
new.pmz <- 0
for (i in 2:length(pmz)) {
  if (pmz[i] == 0) {
    x <-
      0
  } else {
    if (pmz[(i - 1)] == 0) {
      x <-
        peaks(aa, (i - 1))[which.min(abs(pmz[i] - peaks(aa, (i - 1))[, 1])), 1]
    } else {
      if (pmz[(i - 2)] == 0) {
        x <-
          peaks(aa, (i - 2))[which.min(abs(pmz[i] - peaks(aa, (i - 2))[, 1])), 1]
      } else {
        if (pmz[(i - 3)] == 0) {
          x <-
            peaks(aa, (i - 3))[which.min(abs(pmz[i] - peaks(aa, (i - 3))[, 1])), 1]
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

## Create a matrix with precursor m/z and retention time for all spectra in ms2list
precursor <- cbind(new.pmz, header(aa)$retentionTime)
precursor2 <- precursor[ms2log,][vec,]

## Exclude everything with RT >25min
cutend <- precursor2[, 2] < 25 * 60
precursormzrt <- precursor2[cutend,]
ms2list <- ms2spectra2[cutend]

## Get median m/z and median RT for precursor masses that differ less than 10ppm
flist <- list()
mz1 <- precursormzrt
while (nrow(mz1) >= 1) {
  l1 <- abs(mz1[1, 1] - mz1[, 1]) <= mz1[1, 1] * 1E-5
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
### Calculate median m/z and median RT for all 'features'
for (i in 1:length(flist)) {
  if (is.matrix(flist[[i]])) {
    flist[[i]] <- cbind(flist[[i]],
                        rep(median(flist[[i]][, 1]), times = nrow(flist[[i]])),
                        rep(median(flist[[i]][, 2]), times = nrow(flist[[i]])))
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
#### Create IDs
medmzrt$id <- paste("M", round((medmzrt$med.mz), 2),
                    "T", round((medmzrt$med.rt), 2), sep = "")
### Put back in original order
medmzrt <- medmzrt[order(medmzrt$rt),]

## Define merging functions
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

mergeSpecList <- function(speclist, mzmed) {
  mrgls <- list()
  for (z in 1:length(speclist)) {
    z0 <- c()
    for (j in 1:(z - 1)) {
      z0[j] <- mzmed[z] == mzmed[j]
    }
    if (z != 1 & any(z0)) {
      mrgls[[z]] <- mrgls[[which(z0 == T)[1]]]
    } else {
      if (sum(mzmed == mzmed[z]) > 1) {
        z1 <- as.matrix(Reduce(mergeTolerance,
                               speclist[mzmed == mzmed[z]]))
        z1[is.na(z1)] <- 0
        z2 <-
          cbind(z1[, 1], round((rowSums(z1) - z1[, 1]) / ncol(z1[,-1])))
        mrgls[[z]] <- z2
      } else {
        mrgls[[z]] <- speclist[[z]]
      }
    }
  }
  names(mrgls) <- names(speclist)
  mrgls
}

## Merge spectra in the list that fulfill identity criteria
names(ms2list) <- medmzrt$id
mergedlist <- mergeSpecList(ms2list, medmzrt$id)
shortlist <- mergedlist[!duplicated(mergedlist)]
shortmzrt <- medmzrt[!duplicated(mergedlist), 3:4]

## Make list with neutral loss spectra
nllist <- list()
for (i in 1:length(shortlist)) {
  nl <-
    cbind((shortmzrt[i, 1] - shortlist[[i]][, 1]), shortlist[[i]][, 2])
  nl <- subset(nl, nl[, 1] >= -(shortmzrt[i, 1] * 1e-5))
  nllist[[i]] <- nl
}
names(nllist) <- names(shortlist)

## Print precursor m/z and RT from all merged spectra
## and identify in Bruker DataAnalysis
write.table(
  cbind(names(shortlist), shortmzrt),
  file = "161019metaboident_SE_pre.csv",
  sep = ",",
  row.names = F
)

## Read in manual annotations
ident <-
  read.csv(file = "161022metaboident_SE_post.csv", stringsAsFactors = F)
metabonames <- c()
for (n in 1:nrow(ident)) {
  if (is.na(ident[n, 4]) && is.na(ident[n, 5])) {
    metabonames[n] <- ident[n, 1]
  }
  if (!is.na(ident[n, 4])) {
    metabonames[n] <- paste(ident[n, 1], " - ", ident[n, 4], sep = "")
  }
  if (is.na(ident[n, 4]) && !is.na(ident[n, 5])) {
    metabonames[n] <-
      paste(ident[n, 1], " - ", "(", ident[n, 5], ")", sep = "")
  }
}
names(shortlist) <- metabonames
names(nllist) <- metabonames

## Define similarity score
cossim <- function(x, y) {
  colnames(x) <- NULL
  colnames(y) <- NULL
  mm <- mergeTolerance(x, y)
  sum(sqrt(mm[, 2]) * sqrt(mm[, 3])) / (sqrt(sum(mm[, 2])) * sqrt(sum(mm[, 3])))
}

## Create distance matrix for MS2 spectra
distmat <-
  matrix(nrow = length(shortlist), ncol = length(shortlist))
for (m in 1:length(shortlist)) {
  for (l in m:length(shortlist)) {
    if (is.na(distmat[m, l])) {
      distmat[m, l] <- 1 - cossim(shortlist[[m]], shortlist[[l]])
      distmat[l, m] <- distmat[m, l]
    }
  }
}
colnames(distmat) <- names(shortlist)
row.names(distmat) <- names(shortlist)

## Create distance matrix for neutral loss spectra
distmat.nl <-
  matrix(nrow = length(nllist), ncol = length(nllist))
for (m in 1:length(nllist)) {
  for (l in m:length(nllist)) {
    if (is.na(distmat.nl[m, l])) {
      distmat.nl[m, l] <- 1 - cossim(nllist[[m]], nllist[[l]])
      distmat.nl[l, m] <- distmat.nl[m, l]
    }
  }
}
distmat.nl[is.na(distmat.nl)] <- 1
colnames(distmat.nl) <- names(nllist)
row.names(distmat.nl) <- names(nllist)

## Multidimensional scaling

fit <- cmdscale(as.dist(distmat), eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
pdf(file = "figure_mds.pdf",
    height = 12,
    width = 12)
plot(
  x,
  y,
  xlab = "Coordinate 1",
  ylab = "Coordinate 2",
  type = "p",
  col = grey(0.3)
)
dev.off()

## Density based clustering using OPTICS
opt <-
  optics(as.dist(distmat),
         eps = 10000,
         minPts = 3,
         search = "dist")
opt.nl <-
  optics(
    as.dist(distmat.nl),
    eps = 10000,
    minPts = 3,
    search = "dist"
  )

### Identify clusters by cutting the reachability plot (black is noise)
res <- optics_cut(opt, eps_cl = 0.5)

pdf(file = "20161025_optics_SE.pdf",
    height = 6,
    width = 12)
plot(res)
dev.off()

clustmat <- NULL
for (i in c(1:max(res$cluster), 0)) {
  x <-
    cbind(colnames(distmat)[res$cluster == i], rep(i, length(colnames(distmat)[res$cluster == i])))
  clustmat <- rbind(clustmat, x)
}
write.csv(clustmat, file = "20161024_dbclust.csv")
res.nl <- optics_cut(opt.nl, eps_cl = 0.7)

### Create plot
opal <- palette()
palette(c(opal, rep(c("orange", opal[-1]),10)))

pdf(file = "figure_optics.pdf",
    height = 6,
    width = 12)
plot(res)
dev.off()
palette(opal)

clustmat.nl <- NULL
for (i in c(1:max(res.nl$cluster), 0)) {
  x <-
    cbind(colnames(distmat)[res.nl$cluster == i], rep(i, length(colnames(distmat.nl)[res.nl$cluster == i])))
  clustmat.nl <- rbind(clustmat.nl, x)
}
write.csv(clustmat.nl, file = "20161024_dbclustnl.csv")

## Hierarchical clustering
clust <- hclust(as.dist(distmat), method = "average")
hclusttree <- cutree(clust, h = 0.95)
hclustmat <- cbind(names(hclusttree), hclusttree)

### Plot heatmap
hm <- heatmap(distmat, Rowv = as.dendrogram(clust), Colv = "Rowv", distfun = NULL, symm = T)

### Plot dendrogram (2 different layouts)
clr <- brewer.pal(n = 8, name = "Dark2")
pdf(file = "figure_hclust_dendrogram.pdf",
    width = 30,
    height = 30)
plot(
  as.phylo(clust),
  type = "fan",
  cex = 0.7,
  tip.color = rep(clr, 16)[hclusttree]
)
dev.off()
pdf(file = "figure_clusters.pdf",
    width = 10,
    height = 35)
plot(
  as.phylo(clust),
  type = "phylo",
  cex = 0.4,
  tip.color = rep(clr, 16)[hclusttree]
)
dev.off()

write.csv(hclusttree, file = "20161024_hclust.csv")

clust.nl <- hclust(as.dist(distmat.nl), method = "average")
hclusttree.nl <- cutree(clust.nl, h = 0.95)
write.csv(hclusttree.nl, file = "20161024_hclustnl.csv")

### Plot dendrogram
pdf(file = "figure_hclust_dendrogram_nl.pdf",
    width = 30,
    height = 30)
plot(
  as.phylo(clust.nl),
  type = "fan",
  cex = 0.7,
  tip.color = rep(clr, 16)[hclusttree.nl]
)
dev.off()

close(aa)

###################

## Tools for analysis & interpretation

### Print spectrum as table from 'shortlist'
print.matrix <- function(m){
  write.table(format(m, justify="right"),
              row.names=F, col.names=F, quote=F)
}
ms2 <- function(x){print.matrix(shortlist[[x]])}

### Plot spectrum from 'shortlist'
specplot <- function(n, list = shortlist) {
  plot(x = list[[n]][,1],
       y = list[[n]][,2] / max(list[[n]][,2]),
       type = "h",
       xlim = c(0, (max(list[[n]][, 1]) * 1.1)),
       xaxs = "i",
       xlab = expression(italic(m/z)),
       ylim = c(0, 1.1),
       yaxs = "i",
       ylab = "intensity relative to base peak",
       main = names(list[n]))
  text(x = (list[[n]][,1])[(list[[n]][,2] / max(list[[n]][,2])) > 0.1],
       y = (list[[n]][,2] / max(list[[n]][,2]))[(list[[n]][,2] / max(list[[n]][,2])) > 0.1],
       labels = round((list[[n]][, 1])[(list[[n]][,2] / max(list[[n]][,2])) > 0.1], 4),
       pos = 3,
       cex = 0.75)
}

### Create mirror plot of two spectra from 'shortlist'
specplot2 <- function(n, o, list = shortlist) {
  plot(
    x = list[[n]][, 1],
    y = list[[n]][, 2] / max(list[[n]][, 2]),
    type = "h",
    xlim = c(0, (max(c(
      list[[n]][, 1], list[[o]][, 1]
    )) * 1.1)),
    xaxs = "i",
    xlab = expression(italic(m / z)),
    ylim = c(-1.2, 1.2),
    yaxs = "i",
    yaxt = "n",
    ylab = "intensity relative to base peak",
    main = paste(names(list[n]), "--", names(list[o]))
  )
  points(x = list[[o]][, 1],
         y = -(list[[o]][, 2] / max(list[[o]][, 2])),
         type = "h")
  abline(a = 0, b = 0)
  axis(2,
       at = seq(-1, 1, 0.5),
       labels = c(1.0, 0.5, 0.0, 0.5, 1.0))
  text(
    x = (list[[n]][, 1])[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1],
    y = (list[[n]][, 2] / max(list[[n]][, 2]))[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1],
    labels = round((list[[n]][, 1])[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1], 4),
    pos = 3,
    cex = 0.75
  )
  text(
    x = (list[[o]][, 1])[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1],
    y = -((list[[o]][, 2] / max(list[[o]][, 2]))[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1]),
    labels = round((list[[o]][, 1])[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1], 4),
    pos = 1,
    cex = 0.75
  )
}