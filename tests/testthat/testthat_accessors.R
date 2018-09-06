library(CluMSID)
context("Accessor functions for MS2spectrum objects")

load(file = system.file("extdata",
                        "annotatedSpeclist.RData",
                        package = "CluMSID"))

test_that("access_id gives a correct output", {
    expect_equal(access_id(annotatedSpeclist[[1]]), "M146.17T59.35")
})

test_that("access_annotation gives a correct output", {
    expect_equal(access_annotation(annotatedSpeclist[[1]]), "spermidine")
})

test_that("access_precursor gives a correct output", {
    expect_equal(access_precursor(annotatedSpeclist[[1]]), 146.1653)
})

test_that("access_rt gives a correct output", {
    expect_equal(access_rt(annotatedSpeclist[[1]]), 59.35)
})

test_that("access_polarity gives a correct output", {
    expect_equal(access_polarity(annotatedSpeclist[[1]]), "positive")
})

test_that("access_spectrum gives a correct output", {
    expect_is(access_spectrum(annotatedSpeclist[[1]]), "matrix")
})

test_that("access_neutral_losses gives a correct output", {
    expect_is(access_neutral_losses(annotatedSpeclist[[1]]), "matrix")
})
