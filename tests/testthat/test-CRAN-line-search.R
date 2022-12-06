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
