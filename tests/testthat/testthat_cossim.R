library(CluMSID)
context("Cosine similarity function")

spec1 <- matrix(  c(50,  800,
                    100, 500,
                    120, 700,
                    180, 100,
                    200, 1000),
                    ncol = 2,
                    byrow = TRUE)

spec2 <- matrix(  c(50,  800,
                    100, 500,
                    120, 700,
                    180, 100,
                    200, 1000),
                    ncol = 2,
                    byrow = TRUE)

spec3 <- matrix(  c(50,  400,
                    100, 250,
                    120, 350,
                    180, 50,
                    200, 500),
                    ncol = 2,
                    byrow = TRUE)

spec4 <- matrix(  c(55,  800,
                    110, 500,
                    132, 700,
                    198, 100,
                    220, 1000),
                    ncol = 2,
                    byrow = TRUE)

spec5 <- methods::new("MS2spectrum",
                    id = "spec5",
                    spectrum = matrix(  c(50,  800,
                                            100, 500,
                                            120, 700,
                                            180, 100,
                                            200, 1000),
                                        ncol = 2,
                                        byrow = TRUE))

test_that("cosine similarity is one for identical spectra", {
    expect_equal(cossim(spec1, spec2), 1)
})

test_that("cosine similarity is one for identical spectra
            with different intensity scale", {
    expect_equal(cossim(spec2, spec3), 1)
})

test_that("cosine similarity is zero for spectra
            with no shared peaks", {
    expect_equal(cossim(spec2, spec4), 0)
})

test_that("cossim throws an error if spectra have different classes", {
    expect_error(cossim(spec1, spec5))
})
