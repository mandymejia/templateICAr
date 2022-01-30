test_that("Data normalization completes", {
  nr <- 170
  nc <- 1000
  x <- norm_BOLD(matrix(rnorm(nr*nc), nrow=nr))
  testthat::expect_true(all(dim(x) == c(nr, nc)))
})
