library(CluMSID)
context("getSpectrum() and related accessory functions")

load(file = system.file("extdata",
                        "annotatedSpeclist.RData",
                        package = "CluMSID"))

test_that("getSpectrum gives a correct output if
          query matches several spectra", {
    expect_is(getSpectrum(annotatedSpeclist,
                            "precursor", 286.18,
                            mz.tol = 1E-03),
                "list")
    expect_equal(length(getSpectrum(annotatedSpeclist,
                                    "precursor", 286.18,
                                    mz.tol = 1E-03)),
                 4)
    expect_is(getSpectrum(annotatedSpeclist,
                          "precursor", 286.18,
                          mz.tol = 1E-03)[[1]],
              "MS2spectrum")
})

test_that("getSimilarities correctly excludes non-hits",{
    expect_gt(min(getSimilarities(annotatedSpeclist[[137]],
                    annotatedSpeclist, hits_only = TRUE)), 0)
    expect_equal(max(getSimilarities(annotatedSpeclist[[137]],
                    annotatedSpeclist, hits_only = TRUE)), 1)
})


