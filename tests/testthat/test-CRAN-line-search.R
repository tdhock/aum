library(testthat)
intercept.decreasing <- data.frame(
  fp.diff=0, fn.diff=0, intercept=c(0,1,0), slope=0)
test_that("error when intercepts are decreasing", {
  expect_error({
    aum:::aumLineSearch(intercept.decreasing, initialAum=10, maxIterations=10)
  }, "intercepts should be non-decreasing")
})
slope.same <- data.frame(
  fp.diff=0, fn.diff=0, intercept=c(0,1,1), slope=c(0,0,0))
test_that("error when slope same", {
  expect_error({
    aum:::aumLineSearch(slope.same, initialAum=10, maxIterations=10)
  }, "slopes should be increasing for equal intercepts")
})
test_that("line search initial auc correct for tie", {
  bin.diffs <- aum::aum_diffs_binary(c(0,1))
  L <- aum::aum_line_search(bin.diffs, pred.vec=c(0,0))
  expect_equal(L$line_search_result$aucAtStep, 0.5)
  expect_equal(L$line_search_result$aucAfterStep, 1)
})
test_that("line search initial auc correct no tie", {
  bin.diffs <- aum::aum_diffs_binary(c(0,1))
  L <- aum::aum_line_search(bin.diffs, pred.vec=c(2,-2))
  expect_equal(L$line_search_result, rbind(
    data.frame(aum=4, step.size=0, aucAtStep=0, aucAfterStep=0),
    data.frame(aum=0, step.size=2, aucAtStep=0.5, aucAfterStep=1)))
})

three.intersect <- data.frame(
  intercept=c(-1,0,1), 
  slope=c(1, 0, -1),
  fp.diff=c(0.5,0,0.5), 
  fn.diff=c(0,-0.5,-0.5))
aum:::aumLineSearch(three.intersect, initialAum = 2, maxIterations = 3)
