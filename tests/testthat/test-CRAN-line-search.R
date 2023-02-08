library(testthat)
intercept.decreasing <- data.frame(
  fp.diff=0, fn.diff=0, intercept=c(0,1,0), slope=0)
test_that("error when intercepts are decreasing", {
  expect_error({
    aum:::aumLineSearch(intercept.decreasing, maxIterations=10)
  }, "intercepts should be non-decreasing")
})
slope.same <- data.frame(
  fp.diff=0, fn.diff=0, intercept=c(0,1,1), slope=c(0,0,0))
test_that("error when slope same", {
  expect_error({
    aum:::aumLineSearch(slope.same, maxIterations=10)
  }, "slopes should be increasing for equal intercepts")
})
test_that("line search initial auc correct for tie", {
  bin.diffs <- aum::aum_diffs_binary(c(0,1))
  L <- aum::aum_line_search(bin.diffs, pred.vec=c(0,0))
  expected.row <- data.frame(
    aum=0, aum.slope.after=0, step.size=0, auc=0.5, auc.after=1)
  expect_equal(L$line_search_result, expected.row)
})
test_that("line search initial auc correct no tie", {
  bin.diffs <- aum::aum_diffs_binary(c(0,1))
  L <- aum::aum_line_search(bin.diffs, pred.vec=c(2,-2))
  expect_equal(L$line_search_result, rbind(
    data.frame(aum=4, aum.slope.after=-2, step.size=0, auc=0, auc.after=0),
    data.frame(aum=0, aum.slope.after=0, step.size=2, auc=0.5, auc.after=1)))
})

test_that("contrived three way tie computed ok", {
  three.intersect <- data.frame(
    intercept=c(-1,0,1), 
    slope=c(1, 0, -1),
    fp.diff=c(0.5,0,0.5), 
    fn.diff=c(0,-0.5,-0.5))
  L <- aum:::aumLineSearch(three.intersect, maxIterations = 2)
  (expected.df <- rbind(
    data.frame(aum=1, aum.slope.after=-1, step.size=0, auc=3/8, auc.after=3/8),
    data.frame(aum=0, aum.slope.after=0.5, step.size=1, auc=0.5, auc.after=5/8)))
  expect_equal(L, expected.df)
})

test_that("contrived four way tie computed ok", {
  four.intersect <- data.frame(
    intercept=c(-1,1,3,5), 
    slope=c(1,-1,1,-1),
    fp.diff=c(0.5,0,0.5,0), 
    fn.diff=c(0,-0.5,0,-0.5))
  L <- aum:::aumLineSearch(four.intersect, maxIterations = 3)
  (expected.df <- rbind(
    data.frame(aum=3, aum.slope.after=-1, step.size=0, auc=1/4, auc.after=1/4),
    data.frame(aum=2, aum.slope.after=-1, step.size=1, auc=1/2, auc.after=3/4),
    data.frame(aum=0, aum.slope.after=0, step.size=3, auc=7/8, auc.after=1)))
  expect_equal(L, expected.df)
})

test_that("join to three way tie computed ok", {
  four.intersect <- data.frame(
    intercept=c(-3,-1,0,3), 
    slope=c(1,-1,0,-1),
    fp.diff=c(0.5,0,0.5,0), 
    fn.diff=c(0,-0.5,0,-0.5))
  L <- aum:::aumLineSearch(four.intersect, maxIterations = 3)
  expect_equal(L$step.size, c(0,1,3))
})

test_that("join to three way tie then another", {
  four.intersect <- data.frame(
    intercept=c(-3,-1,0,3), 
    slope=c(1,0,0,-1),
    fp.diff=c(0.5,0,0.5,0), 
    fn.diff=c(0,-0.5,0,-0.5))
  L <- aum:::aumLineSearch(four.intersect, maxIterations = 4)
  if(require(ggplot2)){
    ggplot()+
      theme_bw()+
      geom_abline(aes(
        slope=slope, intercept=intercept),
        data=four.intersect)+
      geom_vline(aes(
        xintercept=step.size),
        color="red",
        data=L)+
      coord_cartesian(xlim=c(0, 5), ylim=c(-3,3))
  }
  expect_equal(L$step.size, c(0,2,3,4))
})

test_that("several ties", {
  several.intersect <- data.frame(
    intercept=c(-3,-1,0,3,5,9),
    slope=c(1,0,0,-1,1,-1),
    fp.diff=c(0.5,0,0.25,0,0.25,0), 
    fn.diff=c(0,-0.25,0,-0.25,0,-0.5))
  (L <- aum:::aumLineSearch(several.intersect, maxIterations = 7))
  if(require(ggplot2)){
    ggplot()+
      theme_bw()+
      geom_abline(aes(
        slope=slope, intercept=intercept),
        data=several.intersect)+
      geom_vline(aes(
        xintercept=step.size),
        color="red",
        data=L)+
      coord_cartesian(xlim=c(0, 10), ylim=c(-3,9))
  }
  expect_equal(L$step.size, c(0,2,3,4,6,9,10))
})
