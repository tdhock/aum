library(testthat)
diffs <- function(...){
  L <- list(...)
  for(example.i in seq_along(L)){
    L[[example.i]]$example <- example.i-1L
  }
  do.call(rbind, L)
}
models <- diffs(
  data.frame(fp_diff=c( 1, -1, 1),
             fn_diff=c( 0,  0, 0),
             pred=   c(-1,  0, 1)),
  data.frame(fp_diff=c( 0,  0, 0),
             fn_diff=c( -1, 1, -1),
             pred=   c(-1,  0, 1)))
predictions <- c(0,1)
test_that("noncvx 1fp[-1,0] 1fn[0,1]", {
  L <- aum::aum_interface(models, predictions)
  expect_equal(L$aum, 1)
  expect_equal(L$derivative_mat[,1], c(1,1))
  expect_equal(L$derivative_mat[,2], c(-1,-1))
})

