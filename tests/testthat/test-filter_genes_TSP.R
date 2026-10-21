test_that("Dunn filtering skips undefined features", {
  skip_if_not_installed("dunn.test")

  labels <- rep(c("A", "B", "C"), each = 3)
  x <- rbind(
    varying = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    constant = rep(1, 9),
    with_na = c(1, 2, NA, 4, 5, 6, 7, 8, 9)
  )
  names(x) <- paste0("sample", seq_len(ncol(x)))
  got <- multiclassPairs:::do_dunn_test(x, labels)

  expect_named(got, "varying")
  expect_length(got, 1)
  expect_false(anyNA(got[[1]]$P))
})

test_that("one-vs-one and one-vs-rest filtering handle edge features", {
  skip_if_not_installed("dunn.test")
  skip_if_not_installed("switchBox")

  set.seed(1)
  labels <- rep(c("A", "B", "C"), each = 3)
  x <- matrix(rnorm(27), nrow = 3,
              dimnames = list(c("varying", "constant", "with_na"),
                              paste0("sample", seq_len(9))))
  x["constant", ] <- 1
  x["with_na", 3] <- NA

  object <- multiclassPairs::ReadData(x, labels, verbose = FALSE)

  expect_s3_class(
    multiclassPairs::filter_genes_TSP(
      object, filter = "one_vs_one", featureNo = 2,
      UpDown = FALSE, verbose = FALSE),
    "OnevsrestScheme_genes_TSP"
  )
  expect_s3_class(
    multiclassPairs::filter_genes_TSP(
      object, filter = "one_vs_rest", featureNo = 2,
      UpDown = FALSE, verbose = FALSE),
    "OnevsrestScheme_genes_TSP"
  )
})