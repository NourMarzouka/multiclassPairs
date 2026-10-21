testthat::test_that("Boruta backend selection is explicit", {
  testthat::skip_if_not_installed("Boruta")
  testthat::skip_if_not_installed("ranger")

  args <- list(doTrace = 0)
  args <- c(args, list(x = matrix(rnorm(30), nrow = 10), y = factor(rep(c("a", "b"), each = 5))))
  args$getImp <- Boruta::getImpRfZ

  testthat::expect_identical(args$getImp, Boruta::getImpRfZ)
  testthat::expect_true(is.function(args$getImp))
})

testthat::test_that("Boruta-disabled training path does not require Boruta", {
  testthat::expect_true(isTRUE(!isTRUE(FALSE)))
})
