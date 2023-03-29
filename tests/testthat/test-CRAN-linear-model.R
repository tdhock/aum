library(testthat)
## simulated binary classification problem.
N.rows <- 10
N.cols <- 2
set.seed(1)
feature.mat <- matrix(rnorm(N.rows*N.cols), N.rows, N.cols)
n.folds <- 3
label.vec <- rep(0, N.rows)
label.vec[seq(1, n.folds-1)] <- 1
diffs.dt <- aum::aum_diffs_binary(label.vec)
test_that("error when there are not enough data", {
  expect_error({
    aum::aum_linear_model_cv(feature.mat, diffs.dt, n.folds=n.folds)
  }, "not enough data for 3-fold cross-validation, because there are only 2 examples for which there are non-zero values for the minority diff, fn")
})

label.vec[seq(1, n.folds)] <- 1
diffs.dt <- aum::aum_diffs_binary(label.vec)
test_that("error when there are not enough data", {
  model <- aum::aum_linear_model_cv(feature.mat, diffs.dt, n.folds=n.folds)
  expect_is(model, "aum_linear_model_cv")
})


