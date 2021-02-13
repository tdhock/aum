library(testthat)
diffs <- function(...){
  L <- list(...)
  for(example.i in seq_along(L)){
    L[[example.i]]$example <- example.i-1L
  }
  do.call(rbind, L)
}
models <- diffs(
  data.frame(fp_diff=c( 1,-1, 1),
             fn_diff=c( 0, 0, 0),
             pred=   c(-1, 0, 1)),
  data.frame(fp_diff=c( 0, 0, 0),
             fn_diff=c(-1, 1,-1),
             pred=   c(-1, 0, 1)))
predictions <- c(0,1)
test_that("aum=1 noncvx 1fp[-1,0] 1fn[0,1]", {
  L <- aum::aum_interface(models, predictions)
  expect_equal(L$aum, 1)
  expect_equal(L$derivative_mat[1,], c(-1,0))
  expect_equal(L$derivative_mat[2,], c(-0,1))
})

models <- diffs(
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred=0))
predictions <- c(0,0)
test_that("aum=0 nondiff 1fp[-1,0] 1fn[0,1]", {
  L <- aum::aum_interface(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(-1,0))
  expect_equal(L$derivative_mat[2,], c(0,1))
})

predictions <- c(-1, 1)
test_that("aum=0 diff 1fp[0,0] 1fn[0,0]", {
  L <- aum::aum_interface(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,0))
  expect_equal(L$derivative_mat[2,], c(0,0))
})

predictions <- c(1, -1)
test_that("aum=2 diff 1fp[-1,-1] 1fn[1,1]", {
  L <- aum::aum_interface(models, predictions)
  expect_equal(L$aum, 2)
  expect_equal(L$derivative_mat[1,], c(-1,-1))
  expect_equal(L$derivative_mat[2,], c(1,1))
})
