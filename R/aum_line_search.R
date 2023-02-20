aum_line_search <- structure(function
### Exact line search.
(error.diff.df,
### aum_diffs data frame with B rows, one for each breakpoint in
### example-specific error functions.
  feature.mat,
### N x p matrix of numeric features.
  weight.vec,
### p-vector of numeric linear model coefficients.
  pred.vec=NULL,
### N-vector of numeric predicted values. If NULL, feature.mat and
### weight.vec will be used to compute predicted values.
  maxIterations=nrow(error.diff.df)
### positive int: max number of line search iterations.
){
  . <- fp.diff <- fn.diff <- intercept <- slope <- step.size <- NULL
  ## Above to suppress CRAN NOTE.
  pred.null <- is.null(pred.vec)
  if(pred.null){
    pred.vec <- feature.mat %*% weight.vec
  }
  L <- aum(error.diff.df, pred.vec)
  L$pred.vec <- pred.vec
  L$gradient_pred <- rowMeans(L$derivative_mat)
  L$gradient <- if(pred.null){
    L$gradient_weight <- t(feature.mat) %*% L$gradient_pred
    feature.mat %*% L$gradient_weight
  }else{
    L$gradient_pred
  }
  pred.i <- error.diff.df$example+1L
  L$line_search_input <- data.table(
    fp.diff=error.diff.df$fp_diff,
    fn.diff=error.diff.df$fn_diff,
    intercept=error.diff.df$pred-pred.vec[pred.i],
    slope=L$gradient[pred.i]
  )[, .(
    fp.diff=sum(fp.diff),
    fn.diff=sum(fn.diff)
  ), keyby=.(intercept, slope)]
  line.search.all <- aumLineSearch(L$line_search_input, maxIterations)
  L$line_search_result <- data.table(line.search.all)[0 <= step.size]
  class(L) <- c("aum_line_search", class(L))
  L
### List of class aum_line_search.
}, ex=function(){

  ## Example 1: two binary data.
  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  if(requireNamespace("ggplot2"))plot(bin.diffs)
  bin.line.search <- aum::aum_line_search(bin.diffs, pred.vec=c(10,-10))
  if(requireNamespace("ggplot2"))plot(bin.line.search)

  ## Example 2: two changepoint examples, one with three breakpoints.
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  (nb.diffs <- aum::aum_diffs_penalty(nb.err, c("1.1", "4.2")))
  if(requireNamespace("ggplot2"))plot(nb.diffs)
  nb.line.search <- aum::aum_line_search(nb.diffs, pred.vec=c(1,-1))
  if(requireNamespace("ggplot2"))plot(nb.line.search)
  aum::aum_line_search(nb.diffs, pred.vec=c(1,-1)-c(1,-1)*0.5)

  ## Example 3: all changepoint examples, with linear model.
  X.sc <- scale(neuroblastomaProcessed$feature.mat)
  keep <- apply(is.finite(X.sc), 2, all)
  X.keep <- X.sc[1:50,keep]
  weight.vec <- rep(0, ncol(X.keep))
  (nb.diffs <- aum::aum_diffs_penalty(nb.err, rownames(X.keep)))
  nb.weight.search <- aum::aum_line_search(
    nb.diffs,
    feature.mat=X.keep,
    weight.vec=weight.vec, 
    maxIterations = 200)
  if(requireNamespace("ggplot2"))plot(nb.weight.search)

  ## Alternate viz with x=iteration instead of step size.
  nb.weight.full <- aum::aum_line_search(
    nb.diffs,
    feature.mat=X.keep,
    weight.vec=weight.vec, 
    maxIterations = 1000)
  library(data.table)
  weight.result.tall <- suppressWarnings(melt(
    nb.weight.full$line_search_result[, iteration:=1:.N][, .(
      iteration, auc, q.size,
      log10.step.size=log10(step.size),
      log10.aum=log10(aum))],
    id.vars="iteration"))
  if(require(ggplot2)){
    ggplot()+
      geom_point(aes(
        iteration, value),
        shape=1,
        data=weight.result.tall)+
      facet_grid(variable ~ ., scales="free")+
      scale_y_continuous("")
  }

  ## Example 4: many changepoint examples, optimize predictions.
  all.ids <- rownames(neuroblastomaProcessed$feature.mat)
  all.diffs <- aum::aum_diffs_penalty(nb.err, all.ids)
  current.pred <- rep(0, length(all.ids))
  nb.all.search <- aum::aum_line_search(
    all.diffs, pred.vec=current.pred, maxIterations=2e5)
  some.wide <- data.table(nb.all.search$line_search_result)[
    as.integer(seq(1, .N, l=100))
  ][, log10.aum := log10(aum)]
  some.tall <- melt(some.wide, measure.vars = c("auc", "log10.aum"))
  if(require(ggplot2)){
    ggplot()+
      geom_point(aes(
        step.size, value),
        shape=1,
        data=some.tall)+
      facet_grid(variable ~ ., scales="free")+
      theme_bw()
  }
  
})

plot.aum_line_search <- function
### Plot method for aum_line_search which shows AUM and threshold functions. 
(x,
### list with class "aum_line_search".
  ...
### ignored.
){
  step.size <- aum <- slope <- intercept <- auc.after <- auc <- 
    min.step.size <- max.step.size <- aum.slope.after <- n.data <-
      aum.after <- value <- NULL
  ## Above to suppress CRAN check NOTE.
  aum.df <- data.frame(panel="aum", x$line_search_result)
  last <- x$line_search_result[.N]
  step.after <- last$step*1.05
  step.diff <- step.after-last$step
  last.seg <- data.table(
    panel="aum",
    last,
    step.after, 
    aum.after=last[, step.diff*aum.slope.after+aum])
  auc.segs <- x$line_search_result[, data.table(
    panel="auc", 
    min.step.size=step.size,
    max.step.size=c(step.size[-1],Inf),
    auc=auc.after)]
  all.points <- suppressWarnings(melt(
    x$line_search_result,
    measure.vars=c("auc","intervals","intersections","q.size"),
    variable.name="panel"))
  abline.df <- data.frame(panel="threshold", x$line_search_input)
  hline.df <- data.frame(panel="q.size", n.data=nrow(x$line_search_input))
  ggplot2::ggplot()+
    ggplot2::theme_bw()+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=step.size),
      color="grey",
      data=x$line_search_result)+
    ggplot2::geom_hline(ggplot2::aes(
      yintercept=n.data),
      data=hline.df)+
    ggplot2::geom_text(ggplot2::aes(
      -Inf, n.data, label=sprintf(
        "N=%d", n.data)),
      hjust=0,
      vjust=1.1,
      data=hline.df)+
    ggplot2::geom_line(ggplot2::aes(
      step.size, aum),
      size=1,
      data=aum.df)+
    ggplot2::geom_segment(ggplot2::aes(
      min.step.size, auc,
      xend=max.step.size, yend=auc),
      size=1,
      data=auc.segs)+
    ggplot2::geom_segment(ggplot2::aes(
      step.size, aum,
      xend=step.after, yend=aum.after),
      linetype="dotted",
      size=1,
      data=last.seg)+
    ggplot2::geom_point(ggplot2::aes(
      step.size, value),
      shape=1,
      data=all.points)+
    ggplot2::facet_grid(panel ~ ., scales="free")+
    ggplot2::geom_abline(ggplot2::aes(
      slope=slope, intercept=intercept),
      data=abline.df)+
    ggplot2::geom_point(ggplot2::aes(
      0, intercept),
      data=abline.df)+
    ggplot2::scale_y_continuous("")+
    ggplot2::scale_x_continuous("Step size")
### ggplot.
}

aum_line_search_grid <- structure(function
### Line search for predicted values, with grid search to check.
(error.diff.df,
### aum_diffs data frame with B rows, one for each breakpoint in
### example-specific error functions.
  feature.mat,
### N x p matrix of numeric features.
  weight.vec,
### p-vector of numeric linear model coefficients.
  pred.vec=NULL,
### N-vector of numeric predicted values. If missing, feature.mat and
### weight.vec will be used to compute predicted values.
  maxIterations=nrow(error.diff.df),
### positive int: max number of line search iterations.
  n.grid=10L,
### positive int: number of grid points for checking.
  add.breakpoints=FALSE
### add breakpoints from exact search to grid search.
){
  . <- fp_diff <- fn_diff <- fn_before <- fp_before <-
    FPR_before <- TPR_before <- FPR <- TPR <- NULL
  ## Above for CRAN check NOTE.
  L <- aum_line_search(error.diff.df, feature.mat, weight.vec, pred.vec, maxIterations)
  step.size <- unique(sort(c(
    if(add.breakpoints)L$line_search_result$step.size,
    seq(0, max(L$line_search_result$step.size), l=n.grid))))
  step.mat <- matrix(step.size, length(L$pred.vec), length(step.size), byrow=TRUE)
  pred.mat <- as.numeric(L$pred.vec)-step.mat*as.numeric(L$gradient)
  totals <- colSums(error.diff.df[, .(fp_diff, fn_diff)])
  grid.dt.list <- list()
  for(pred.col in 1:ncol(pred.mat)){
    pred <- pred.mat[,pred.col]
    grid.aum <- aum::aum(error.diff.df, pred)
    before.dt <- data.table(grid.aum$total_error, key="thresh")[, `:=`(
      TPR_before=1-fn_before/-totals[["fn_diff"]],
      FPR_before=fp_before/totals[["fp_diff"]])]
    auc <- before.dt[, .(
      FPR=c(FPR_before, 1),
      TPR=c(TPR_before, 1)
    )][, sum((FPR[-1]-FPR[-.N])*(TPR[-1]+TPR[-.N])/2)]
    grid.dt.list[[pred.col]] <- data.table(
      step.size=step.size[pred.col],
      aum=grid.aum$aum,
      auc)
  }
  L$grid_aum <- rbindlist(grid.dt.list)
  class(L) <- c("aum_line_search_grid", class(L))
  L
### List of class aum_line_search_grid.
}, ex=function(){

  ## Example 1: two binary data.
  (bin.diffs <- aum::aum_diffs_binary(c(1,0)))
  if(requireNamespace("ggplot2"))plot(bin.diffs)
  bin.line.search <- aum::aum_line_search_grid(bin.diffs, pred.vec=c(-10,10))
  if(requireNamespace("ggplot2"))plot(bin.line.search)

  ## Example 2: two changepoint examples, one with three breakpoints.
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  (nb.diffs <- aum::aum_diffs_penalty(nb.err, c("4.2", "1.1")))
  if(requireNamespace("ggplot2"))plot(nb.diffs)
  (nb.line.search <- aum::aum_line_search_grid(nb.diffs, pred.vec=c(-1,1)))
  if(requireNamespace("ggplot2"))plot(nb.line.search)

  ## Example 3: 50 changepoint examples, with linear model.
  X.sc <- scale(neuroblastomaProcessed$feature.mat[1:50,])
  keep <- apply(is.finite(X.sc), 2, all)
  X.keep <- X.sc[,keep]
  weight.vec <- rep(0, ncol(X.keep))
  nb.diffs <- aum::aum_diffs_penalty(nb.err, rownames(X.keep))
  nb.weight.search <- aum::aum_line_search_grid(
    nb.diffs,
    feature.mat=X.keep,
    weight.vec=weight.vec,
    maxIterations = 200)
  if(requireNamespace("ggplot2"))plot(nb.weight.search)

  ## Example 4: counting intersections and intervals at each
  ## iteration/step size, when there are ties.
  (bin.diffs <- aum::aum_diffs_binary(c(0,0,0,1,1,1)))
  bin.line.search <- aum::aum_line_search_grid(
    bin.diffs, pred.vec=c(2,3,-1,1,-2,0), n.grid=21) 
  if(require("ggplot2")){
    plot(bin.line.search)+
      geom_text(aes(
        step.size, Inf, label=sprintf(
          "%d,%d", intersections, intervals)),
        vjust=1.1,
        data=data.frame(
          panel="threshold", bin.line.search$line_search_result))
  }

})
  
plot.aum_line_search_grid <- function
### Plot method for aum_line_search_grid which shows AUM and threshold
### functions, along with grid points for checking.
(x,
### list with class "aum_line_search_grid".
  ...
### ignored.
){
  auc.after <- min.step.size <- max.step.size <- auc <- value <-
    step.size <- aum <- slope <- intercept <- search <- NULL
  ## Above to suppress CRAN check NOTE.
  aum.df <- data.frame(
    search="exact", panel="aum", x$line_search_result)
  abline.df <- data.frame(
    search="exact", panel="threshold", x$line_search_input)
  grid.df <- melt(
    data.table(search="grid", x$grid_aum), 
    measure.vars = c("auc", "aum"), 
    variable.name="panel")
  auc.segs <- x$line_search_result[, data.table(
    panel="auc", 
    min.step.size=step.size,
    max.step.size=c(step.size[-1],Inf),
    auc=auc.after)]
  auc.points <- x$line_search_result[, data.table(
    panel="auc",
    step.size, 
    auc)]
  ggplot2::ggplot()+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.spacing=grid::unit(0,"lines"))+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=step.size),
      color="grey",
      data=x$line_search_result)+
    ggplot2::geom_segment(ggplot2::aes(
      min.step.size, auc,
      xend=max.step.size, yend=auc),
      data=auc.segs)+
    ggplot2::geom_point(ggplot2::aes(
      step.size, auc),
      data=auc.points)+
    ggplot2::geom_line(ggplot2::aes(
      step.size, aum, color=search),
      data=aum.df)+
    ggplot2::scale_color_manual(
      values=c(exact="black", grid="red"))+
    ggplot2::facet_grid(panel ~ ., scales="free")+
    ggplot2::geom_abline(ggplot2::aes(
      slope=slope, intercept=intercept, color=search),
      data=abline.df)+
    ggplot2::geom_point(ggplot2::aes(
      0, intercept, color=search),
      shape=1,
      data=abline.df)+
    ggplot2::scale_y_continuous("")+
    ggplot2::geom_point(
      ggplot2::aes(
        step.size, value, color=search),
      shape=1,
      data=grid.df)
### ggplot.
}
