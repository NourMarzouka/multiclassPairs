testthat::test_that("train_RF supports Boruta enabled and disabled", {
  testthat::skip_if_not_installed("Boruta")
  testthat::skip_if_not_installed("ranger")

  set.seed(42)
  labels <- factor(rep(c("A", "B", "C"), each = 8))
  x <- matrix(rnorm(6 * length(labels)), nrow = 6,
              dimnames = list(paste0("g", 1:6), NULL))
  x[1, labels == "A"] <- x[1, labels == "A"] + 3
  x[2, labels == "B"] <- x[2, labels == "B"] + 3
  x[3, labels == "C"] <- x[3, labels == "C"] + 3
  colnames(x) <- paste0("s", seq_len(ncol(x)))

  data_object <- ReadData(x, labels, verbose = FALSE)

  rules <- list(
    all = rbind(c("g1", "g2"), c("g3", "g4"), c("g5", "g6")),
    A = rbind(c("g1", "g2"), c("g1", "g3")),
    B = rbind(c("g3", "g4"), c("g2", "g4")),
    C = rbind(c("g5", "g6"), c("g2", "g5"))
  )
  sorted_rules <- list(list(sorted_rules = rules))
  class(sorted_rules) <- "RandomForest_sorted_rules"

  common <- list(
    data_object = data_object,
    sorted_rules_RF = sorted_rules,
    rules_altogether = 0,
    rules_one_vs_rest = 2,
    num.trees = 10,
    write.forest = TRUE,
    keep.inbag = TRUE,
    probability = TRUE,
    verbose = FALSE
  )

  no_boruta <- do.call(train_RF, c(common, list(run_boruta = FALSE)))
  testthat::expect_s3_class(no_boruta, "rule_based_RandomForest")
  testthat::expect_null(no_boruta$RF_scheme$boruta)

  with_boruta <- do.call(train_RF, c(
    common,
    list(run_boruta = TRUE,
         boruta_args = list(doTrace = 0, maxRuns = 5))
  ))
  testthat::expect_s3_class(with_boruta, "rule_based_RandomForest")
  testthat::expect_s3_class(with_boruta$RF_scheme$boruta, "Boruta")
  testthat::expect_identical(with_boruta$RF_scheme$boruta$call$getImp,
                             Boruta::getImpRfZ)
})
