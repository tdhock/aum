library(testthat)
library(data.table)

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
  example=0:1,
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
  example=0:1,
  pred=0,
  fp_diff=c(1, 0),
  fn_diff=c(0, -1),
  row.names=NULL,
  stringsAsFactors=FALSE)
test_that("binary diffs computed for three named labels", {
  (computed <- aum::aum_diffs_binary(c(a=0,b=1,c=0), c("c","b")))
  expect_equal(as.data.frame(computed[order(example)]), exp.df)
})

test_that("error for numeric example", {
  expect_error({
    aum::aum_diffs(1, 1, 1, 1)
  }, "example must be integer vector but has class: numeric")
})

test_that("error for non-unique predicted example names", {
  expect_error({
    aum::aum_diffs("ex1", 1, 1, 1, c("a","a","b","b","c"))
  }, "elements of pred.name.vec must be unique, problems: a, b")
})

test_that("error for columns in penalty error", {
  expect_error({
    aum::aum_diffs_penalty(data.frame(example=1L, min.lambda=0, fp=0))
  }, "errors.df must have numeric column named fn")
  expect_error({
    aum::aum_diffs_penalty(data.frame(example=1L, min.lambda=0, fn=0))
  }, "errors.df must have numeric column named fp")
  expect_error({
    aum::aum_diffs_penalty(data.frame(example=1L, fp=0, fn=0))
  }, "errors.df must have numeric column named min.lambda")
  expect_error({
    aum::aum_diffs_penalty(data.frame(example=1, min.lambda=0, fp=0, fn=0))
  }, "errors.df must have integer or character column named example")
})

test_that("error if min.lambda does not start at 0", {
  simple.df <- data.frame(
    example=1L,
    min.lambda=exp(1:4),
    fp=c(10,4,4,0),
    fn=c(0,2,2,10))
  expect_error({
    aum::aum_diffs_penalty(simple.df, denominator="count")
  }, "need min.lambda=0 for each example, problems: 1")
})

test_that("error if min.lambda repeated", {
  simple.df <- data.frame(
    example=1L,
    min.lambda=c(0, 1, 1, 2),
    fp=c(10,4,4,0),
    fn=c(0,2,2,10))
  expect_error({
    aum::aum_diffs_penalty(simple.df, denominator="count")
  },
  "need only one min.lambda per example, problems with more are (example:min.lambda) 1:1",
  fixed=TRUE)
})

test_that("rate works for one ex", {
  simple.df <- data.frame(
    example=1L,
    min.lambda=c(0, exp(1:3)),
    fp=c(10,4,4,0),
    fn=c(0,2,2,10))
  (simple.diffs <- aum::aum_diffs_penalty(simple.df, denominator="count"))
  expect_equal(simple.diffs$pred, c(-3, -1))
  expect_equal(simple.diffs$fp_diff, c(4, 6))
  expect_equal(simple.diffs$fn_diff, c(-8, -2))
  (simple.rates <- aum::aum_diffs_penalty(simple.df, denominator="rate"))
  expect_equal(simple.rates$pred, c(-3, -1))
  expect_equal(simple.rates$fp_diff, c(0.4, 0.6))
  expect_equal(simple.rates$fn_diff, c(-0.8, -0.2))
})

test_that("rate works for three ex, one with no diffs", {
  four.dt <- rbind(
    data.table(
      example="one",
      min.lambda=c(0, exp(1:3)),
      fp=c(10,4,4,0),
      fn=c(0,2,2,8)),
    data.table(
      example="two",
      min.lambda=c(0, exp(4:6)),
      fp=c(1,0,0,0),
      fn=c(0,0,0,2)),
    data.table(
      example="three",
      min.lambda=c(0, exp(44:46)),
      fp=c(100,0,0,0),
      fn=c(0,0,0,100)),
    data.table(
      example="constantFN",
      min.lambda=c(0,exp(9)),
      fp=c(11,2),
      fn=c(1,1)))
  three.ids <- c("one","constantFN","two")
  (count.diffs <- aum::aum_diffs_penalty(four.dt, three.ids, denominator="count"))
  expected.counts <- data.frame(
    example=c(0,0,1,2,2),
    pred=c(-3,-1,-9,-6,-4),
    fp_diff=c(4,6,9,0,1),
    fn_diff=c(-6,-2,0,-2,0))
  expect_equal(data.frame(count.diffs), expected.counts)
  (rate.diffs <- aum::aum_diffs_penalty(four.dt, three.ids, denominator="rate"))
  (expected.rates <- with(expected.counts, data.frame(
    example, pred, fp_diff=fp_diff/sum(fp_diff), fn_diff=-fn_diff/sum(fn_diff))))
  expect_equal(data.frame(rate.diffs), expected.rates)
})

test_that("aum_errors works even if input not sorted", {
  diff.df <- data.frame(
    example=as.integer(c(0, 0, 1)),
    pred=c(2, 1, 3),
    fp_diff=c(1, 0, 1),
    fn_diff=c(0, -1, 0))
  computed.dt <- aum::aum_errors(diff.df)
  exp.dt <- data.table(
    example=as.integer(c(0,0,0,1,1)),
    min.pred=c(-Inf,1,2,-Inf,3),
    max.pred=c(1,2,Inf,3,Inf),
    fp=c(0,0,1,0,1),
    fn=c(1,0,0,0,0))
  expect_equal(computed.dt, exp.dt)
})
