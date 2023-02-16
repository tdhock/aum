library(testthat)
library(data.table)
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

test_that("contrived three way tie computed ok", {
  three.intersect <- data.frame(
    intercept=c(-1,0,1), 
    slope=c(1, 0, -1),
    fp.diff=c(0.5,0,0.5), 
    fn.diff=c(0,-0.5,-0.5))
  L <- aum:::aumLineSearch(three.intersect, maxIterations = 2)
  (expected.df <- rbind(
    data.frame(
      step.size=0, aum=1, aum.slope.after=-1, 
      auc=3/8, auc.after=3/8,
      intersections=0, intervals=0, q.size=1),
    data.frame(
      step.size=1, aum=0, aum.slope.after=0.5, 
      auc=0.5, auc.after=5/8,
      intersections=1, intervals=2, q.size=0)))
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
    data.frame(
      step.size=0, aum=3, aum.slope.after=-1, 
      auc=1/4, auc.after=1/4,
      intersections=0, intervals=0, q.size=1),
    data.frame(
      step.size=1, aum=2, aum.slope.after=-1, 
      auc=1/2, auc.after=3/4,
      intersections=2, intervals=2, q.size=1),
    data.frame(
      step.size=3, aum=0, aum.slope.after=0, 
      auc=7/8, auc.after=1,
      intersections=1, intervals=1, q.size=0)))
  expect_equal(L, expected.df)
})

test_that("join to three way tie computed ok", {
  four.intersect <- data.frame(
    intercept=c(-3,-1,0,3), 
    slope=c(1,-1,0,-1),
    fp.diff=c(0.5,0,0.5,0), 
    fn.diff=c(0,-0.5,0,-0.5))
  expected.step <- c(0,1,3)
  L <- aum:::aumLineSearch(
    four.intersect, maxIterations = length(expected.step))
  expect_equal(L$step.size, expected.step)
})

test_that("join to three way tie then another", {
  four.intersect <- data.frame(
    intercept=c(-3,-1,0,3), 
    slope=c(1,0,0,-1),
    fp.diff=c(0.5,0,0.5,0), 
    fn.diff=c(0,-0.5,0,-0.5))
  expected.step <- c(0,2,3,4)
  L <- aum:::aumLineSearch(
    four.intersect, maxIterations = length(expected.step))
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
  expect_equal(L$step.size, expected.step)
})

test_that("several ties", {
  several.intersect <- data.frame(
    intercept=c(-3,-1,0,3,5,9),
    slope=c(1,0,0,-1,1,-1),
    fp.diff=c(0.5,0,0.25,0,0.25,0), 
    fn.diff=c(0,-0.25,0,-0.25,0,-0.5))
  expected.step <- c(0,2,3,4,6,9,10)
  (L <- aum:::aumLineSearch(
    several.intersect, maxIterations = length(expected.step)))
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
  expect_equal(L$step.size, expected.step)
})

test_that("2 binary line search ok windows", {
  bin.diffs <- aum::aum_diffs_binary(c(1,0))
  bin.line.search <- aum::aum_line_search(bin.diffs, pred.vec=c(-10,10))
  expected.dt <- rbind(
    data.table(
      step.size=0, aum=20, aum.slope.after=-2, 
      auc=0, auc.after=0,
      intersections=0, intervals=0, q.size=1),
    data.table(
      step.size=10, aum=0, aum.slope.after=0, 
      auc=0.5, auc.after=1,
      intersections=1, intervals=1, q.size=0))
  expect_equal(bin.line.search$line_search_result, expected.dt)
})

test_that("line search initial auc correct for tie", {
  bin.diffs <- aum::aum_diffs_binary(c(0,1))
  L <- aum::aum_line_search(bin.diffs, pred.vec=c(0,0))
  expected.row <- data.table(
    step.size=0, aum=0, aum.slope.after=0, auc=0.5, auc.after=1, 
    intersections=0, intervals=0, q.size=0)
  expect_equal(L$line_search_result, expected.row)
})

test_that("complex real data example", {
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  all.ids <- rownames(neuroblastomaProcessed$feature.mat)
  all.diffs <- aum::aum_diffs_penalty(nb.err, all.ids)
  current.pred <- rep(0, length(all.ids))
  nb.search <- aum::aum_line_search_grid(
    all.diffs, pred.vec=current.pred, maxIterations=2e5)
  expect_true(all(nb.search$line_search_result$aum >= 0))
  some=data.table(
    nb.search$line_search_input
  )[, id := seq(0,.N-1)][J(id=560:565), on="id"]
  points.dt <- some[some, .(
    x.id, i.id, x.intercept, x.slope, i.intercept, i.slope
  ), on=.(id > id), nomatch=0L
  ][, step := (x.intercept-i.intercept)/(i.slope-x.slope)
    ][, thresh := step*x.slope+x.intercept
      ][is.finite(step) & step>0][order(step)]
  if(require(ggplot2)){
    plot(nb.search)
    ggplot()+
      geom_abline(aes(
        slope=slope, intercept=intercept),
        data=some)+
      geom_label(aes(
        0, intercept, label=id),
        data=some)+
      geom_point(aes(
        step, thresh),
        data=points.dt)+
      geom_vline(aes(
        xintercept=step),
        data=data.frame(step=c(0.010611,0.010801)))
  }
  LDF <- aum:::aumLineSearch(some, nrow(points.dt)+1)
  step.dt <- data.table(
    computed=LDF$step, 
    expected=c(0, points.dt$step),
    x.id=c(NA, points.dt$x.id),
    i.id=c(NA, points.dt$i.id))
  step.dt[, expect_equal(computed, expected)]
})

