library(CluMSID)
context("Accessor functions for MS2spectrum objects")

load(file = system.file("extdata",
                        "annotatedSpeclist.RData",
                        package = "CluMSIDdata"))

test_that("accessID gives a correct output", {
    expect_equal(accessID(annotatedSpeclist[[1]]), "M146.17T59.35")
})

test_that("accessAnnotation gives a correct output", {
    expect_equal(accessAnnotation(annotatedSpeclist[[1]]), "spermidine")
})

test_that("accessPrecursor gives a correct output", {
    expect_equal(accessPrecursor(annotatedSpeclist[[1]]), 146.1653)
})

test_that("accessRT gives a correct output", {
    expect_equal(accessRT(annotatedSpeclist[[1]]), 59.35)
})

test_that("accessPolarity gives a correct output", {
    expect_equal(accessPolarity(annotatedSpeclist[[1]]), "positive")
})

test_that("accessSpectrum gives a correct output", {
    expect_is(accessSpectrum(annotatedSpeclist[[1]]), "matrix")
})

test_that("accessNeutralLosses gives a correct output", {
    expect_is(accessNeutralLosses(annotatedSpeclist[[1]]), "matrix")
})
