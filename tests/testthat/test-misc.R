test_that("Data normalization completes", {
  nv <- 1000
  nt <- 170
  x <- norm_BOLD(matrix(rnorm(nv*nt), nrow=nv))
  testthat::expect_true(all(dim(x) == c(nv, nt)))
})
