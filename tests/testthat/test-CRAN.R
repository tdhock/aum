library(testthat)
ggcost <- function(){
  if(interactive() && require("ggplot2") && require("data.table")){
    prob.vec <- unique(models$example)
    gg.dt <- data.table(problem=prob.vec+1L)[, {
      data.table(log.pen=seq(-2, 2, by=0.5))[, {
        predictions[problem] <- log.pen
        L <- aum::aum(models, predictions)
        with(L, data.table(aum))
      }, by=log.pen]
    }, by=problem]
    ggplot()+
      geom_vline(aes(
        xintercept=predictions),
        data=data.table(predictions, problem=seq_along(predictions)))+
      geom_point(aes(
        log.pen, aum),
        data=gg.dt)+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(problem ~ .)+
      coord_equal()
  }
}
diffs <- function(...){
  L <- list(...)
  for(example.i in seq_along(L)){
    L[[example.i]]$example <- example.i-1L
  }
  do.call(rbind, L)
}


models <- diffs(
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred=0))
predictions <- c(a=1, b=-1)
ggcost()
test_that("aum=2 diff 1fp[1,1] 1fn[-1,-1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 2)
  expect_identical(rownames(L$derivative_mat), c("a","b"))
  expect_equal(L$derivative_mat[1,], c(1,1))
  expect_equal(L$derivative_mat[2,], c(-1,-1))
})

predictions <- c(0,0)
ggcost()
test_that("aum=0 nondiff 1fp[0,1] 1fn[-1,0]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,1))
  expect_equal(L$derivative_mat[2,], c(-1,0))
})

predictions <- c(-1, 1)
ggcost()
test_that("aum=0 diff 1fp[0,0] 1fn[0,0]", {
  (L <- aum::aum(models, predictions))
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,0))
  expect_equal(L$derivative_mat[2,], c(0,0))
})

models <- diffs(
  data.frame(fp_diff=c( 1,-1, 1),
             fn_diff=c( 0, 0, 0),
             pred=   c(-1, 0, 1)),
  data.frame(fp_diff=c( 0, 0, 0),
             fn_diff=c(-1, 1,-1),
             pred=   c(-1, 0, 1)))
predictions <- c(0,1)
ggcost()
test_that("aum=1 noncvx 1fp[1,-1] 1fn[1,-1]", {
  (L <- aum::aum(models, predictions))
  expect_equal(L$aum, 1)
  expect_equal(L$derivative_mat[1,], c(1, -1))
  expect_equal(L$derivative_mat[2,], c(1, -1))
})

predictions <- c(0,0)
ggcost()
test_that("aum=0 1fp[-1,2] 1fn[-2,1]", {
  (L <- aum::aum(models, predictions))
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(-1, 2))
  expect_equal(L$derivative_mat[2,], c(-2, 1))
})

predictions <- c(-1,1)
ggcost()
test_that("aum=0 1fp[0,1] 1fn[-1,0]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0, 1))
  expect_equal(L$derivative_mat[2,], c(-1, 0))
})

models <- diffs(
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred   =0),
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0))
predictions <- c(0,0,0)
ggcost()
test_that("1fp[0,1] 1fn[-1,0] 1fp[0,1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,1))
  expect_equal(L$derivative_mat[2,], c(-1,0))
  expect_equal(L$derivative_mat[3,], c(0,1))
})

models <- diffs(
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-2,
             pred   =0),
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =0))
predictions <- c(0,0,0)
ggcost()
test_that("1fp[0,1] 2fn[-2,0] 1fp[0,1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,1))
  expect_equal(L$derivative_mat[2,], c(-2,0))
  expect_equal(L$derivative_mat[3,], c(0,1))
})

models <- diffs(
  data.frame(fp_diff=2,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred   =0),
  data.frame(fp_diff=c(0,1,1),
             fn_diff=c(-1,-1,0),
             pred   =c(-1,0,1)))
predictions <- c(0,0,0)
ggcost()
test_that("2fp[0,2] 1fn[-1,0] 2fp2fn(0)[-1,1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,2))
  expect_equal(L$derivative_mat[2,], c(-1,0))
  expect_equal(L$derivative_mat[3,], c(-1,1))
})

predictions <- c(0,0,-1)
ggcost()
test_that("2fp[1,2] 1fn[-1,0] 2fp2fn(1)[-2,-1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 1)
  expect_equal(L$derivative_mat[1,], c(1,2))
  expect_equal(L$derivative_mat[2,], c(-1,0))
  expect_equal(L$derivative_mat[3,], c(-2,-1))
})

models <- diffs(
  data.frame(fp_diff=4,
             fn_diff=0,
             pred   =0),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred   =0),
  data.frame(fp_diff=c(0,1,1),
             fn_diff=c(-1,-1,0),
             pred   =c(-1,0,1)))
predictions <- c(1,0,0)
ggcost()
test_that("4fp[2,3] 1fn[-1,-1] 2fp2fn[-2,-1]", {
  L <- aum::aum(models, predictions)
  expect_equal(L$derivative_mat[1,], c(2,3))
  expect_equal(L$derivative_mat[2,], c(-1,-1))
  expect_equal(L$derivative_mat[3,], c(-2,-1))
})

models <- diffs(
  data.frame(fp_diff=0,
             fn_diff=1,
             pred   =0))
predictions <- c(1,0,0)
test_that("error for fn<0", {
  expect_error({
    aum::aum(models, predictions)
  }, "fn should be non-negative")
})

models <- diffs(
  data.frame(fp_diff=-1,
             fn_diff=0,
             pred   =0))
predictions <- c(1,0,0)
test_that("error for fp<0", {
  expect_error({
    aum::aum(models, predictions)
  }, "fp should be non-negative")
})

models <- diffs(
  data.frame(fp_diff=1,
             fn_diff=0,
             pred   =-1),
  data.frame(fp_diff=0,
             fn_diff=-1,
             pred=1))
predictions <- c(0,0,0)
test_that("extra pred ok", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 2)
  expect_equal(L$derivative_mat[1,], c(1,1))
  expect_equal(L$derivative_mat[2,], c(-1,-1))
  expect_equal(L$derivative_mat[3,], c(0,0))
  other <- with(models, data.frame(example=example+1L, fp_diff, fn_diff, pred))
  L <- aum::aum(other, predictions)
  expect_equal(L$aum, 2)
  expect_equal(L$derivative_mat[1,], c(0,0))
  expect_equal(L$derivative_mat[2,], c(1,1))
  expect_equal(L$derivative_mat[3,], c(-1,-1))
})

test_that("error for example=length(predictions)", {
  other <- with(models, data.frame(example=example+2L, fp_diff, fn_diff, pred))
  expect_error({
    aum::aum(other, predictions)
  }, "example should be less than number of predictions")
})

test_that("error for example<0", {
  other <- with(models, data.frame(example=example-1L, fp_diff, fn_diff, pred))
  expect_error({
    aum::aum(other, predictions)
  }, "example should be non-negative")
})

test_that("error for no predictions", {
  expect_error({
    aum::aum(models, numeric())
  }, "need at least one prediction")
})

test_that("error for non-finite prediction", {
  expect_error({
    aum::aum(models, c(0,NA))
  }, "all predictions should be finite")
  expect_error({
    aum::aum(models, c(0,Inf))
  }, "all predictions should be finite")
})

models <- data.frame(
  example=integer(),
  fp_diff=numeric(),
  fn_diff=numeric(),
  pred   =numeric())
predictions <- c(1,0,0)
test_that("no diffs ok", {
  L <- aum::aum(models, predictions)
  expect_equal(L$aum, 0)
  expect_equal(L$derivative_mat[1,], c(0,0))
  expect_equal(L$derivative_mat[2,], c(0,0))
  expect_equal(L$derivative_mat[3,], c(0,0))
})

data(neg.zero.fp, package="aum")
test_that("no error for negative zero fp", {
  result <- with(neg.zero.fp, aum::aum(diffs, pred))
  expect_is(result, "list")
})
