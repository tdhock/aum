library(testthat)

test_that("duplicate names is an error", {
  expect_error({
    aum::aum_diffs_binary(c(a=0,b=1,b=0))
  }, "if label.vec has names they must be unique, problems: b")
})

test_that("non-numeric labels is an error", {
  expect_error({
    aum::aum_diffs_binary(factor(c(-1,1)))
  }, "label.vec must be numeric vector with length>0 and all elements either 0,1 or -1,1")
})

test_that("non-finite label is an error", {
  expect_error({
    aum::aum_diffs_binary(c(-1,1,NA))
  }, "label.vec must be numeric vector with length>0 and all elements either 0,1 or -1,1")
})

test_that("non-finite label is an error", {
  expect_error({
    aum::aum_diffs_binary(numeric())
  }, "label.vec must be numeric vector with length>0 and all elements either 0,1 or -1,1")
})

test_that("non-binary label is an error", {
  expect_error({
    aum::aum_diffs_binary(c(-1,1,2))
  }, "label.vec must be numeric vector with length>0 and all elements either 0,1 or -1,1")
})

exp.df <- data.frame(
  example=1:2,
  pred=0,
  fp_diff=c(1, 0),
  fn_diff=c(0, -1),
  row.names=NULL)
test_that("binary diffs computed for two un-named labels", {
  (computed <- aum::aum_diffs_binary(c(0,1)))
  expect_equal(as.data.frame(computed), exp.df)
  (computed <- aum::aum_diffs_binary(c(-1,1)))
  expect_equal(as.data.frame(computed), exp.df)
})

exp.df <- data.frame(
  example=c("a","b","c"),
  pred=0,
  fp_diff=c(1, 0, 1),
  fn_diff=c(0, -1, 0),
  row.names=NULL,
  stringsAsFactors=FALSE)
test_that("binary diffs computed for three named labels", {
  (computed <- aum::aum_diffs_binary(c(a=0,b=1,c=0)))
  expect_equal(as.data.frame(computed), exp.df)
})

