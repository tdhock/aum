aum_line_search <- structure(function
### Line search for predicted values.
(error.diff.df,
### aum_diffs data frame.
  pred.vec,
### Numeric vector of predicted values.
  maxIterations=nrow(error.diff.df)
### positive int: max number of line search iterations.
){
  L <- aum(error.diff.df, pred.vec)
  L$gradient <- rowMeans(L$derivative_mat)
  pred.i <- error.diff.df$example+1L
  L$line_search_input <- data.frame(
    fp.diff=error.diff.df$fp_diff,
    fn.diff=error.diff.df$fn_diff,
    intercept=error.diff.df$pred-pred.vec[pred.i],
    slope=L$gradient[pred.i])
  line.search.all <- aumLineSearch(L$line_search_input, L$aum, maxIterations)
  keep <- 0 <= line.search.all$step.size 
  L$line_search_result <- line.search.all[keep,]
  class(L) <- c("aum_line_search", class(L))
  L
### List of class aum_line_search.
}, ex=function(){

  ## Example 1: two binary data.
  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  if(requireNamespace("ggplot2"))plot(bin.diffs)
  (bin.line.search <- aum::aum_line_search(bin.diffs, c(10,-10)))
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
  (nb.line.search <- aum::aum_line_search(nb.diffs, c(1,-1)))
  if(requireNamespace("ggplot2"))plot(nb.line.search)

})

plot.aum_line_search <- function
### Plot method for aum_line_search which shows AUM and threshold functions. 
(x,
### list with class "aum_line_search".
  ...
### ignored.
){
  step.size <- aum <- slope <- intercept <- NULL
  ## Above to suppress CRAN check NOTE.
  aum.df <- data.frame(panel="aum", x$line_search_result)
  abline.df <- data.frame(panel="threshold", x$line_search_input)
  ggplot2::ggplot()+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.spacing=grid::unit(0,"lines"))+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=step.size),
      color="grey",
      data=x$line_search_result)+
    ggplot2::geom_point(ggplot2::aes(
      step.size, aum),
      data=aum.df)+
    ggplot2::geom_line(ggplot2::aes(
      step.size, aum),
      size=1,
      data=aum.df)+
    ggplot2::facet_grid(panel ~ ., scales="free")+
    ggplot2::geom_abline(ggplot2::aes(
      slope=slope, intercept=intercept),
      data=abline.df)+
    ggplot2::geom_point(ggplot2::aes(
      0, intercept),
      data=abline.df)+
    ggplot2::scale_y_continuous("")
}

aum_line_search_grid <- structure(function(error.diff.df, pred.vec, maxIterations=nrow(error.diff.df), n.grid=10L){
  L <- aum_line_search(error.diff.df, pred.vec, maxIterations)
  step.size <- seq(0, max(L$line_search_result$step.size), l=n.grid)
  step.mat <- matrix(step.size, length(pred.vec), length(step.size), byrow=TRUE)
  pred.mat <- pred.vec-step.mat*L$gradient
  L$grid_aum <- data.table(
    step.size,
    aum=apply(pred.mat, 2, function(pred)aum::aum(error.diff.df, pred)$aum))
  class(L) <- c("aum_line_search_grid", class(L))
  L
}, ex=function(){

  ## Example 1: two binary data.
  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  if(requireNamespace("ggplot2"))plot(bin.diffs)
  (bin.line.search <- aum::aum_line_search_grid(bin.diffs, c(10,-10)))
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
  (nb.line.search <- aum::aum_line_search_grid(nb.diffs, c(1,-1)))
  if(requireNamespace("ggplot2"))plot(nb.line.search)
  
})
  
plot.aum_line_search_grid <- function(x, ...){
  step.size <- aum <- slope <- intercept <- search <- NULL
  ## Above to suppress CRAN check NOTE.
  aum.df <- data.frame(
    search="exact", panel="aum", x$line_search_result)
  abline.df <- data.frame(
    search="exact", panel="threshold", x$line_search_input)
  grid.df <- data.frame(
    search="grid", panel="aum", x$grid_aum)
  ggplot2::ggplot()+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.spacing=grid::unit(0,"lines"))+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=step.size),
      color="grey",
      data=x$line_search_result)+
    ggplot2::geom_point(ggplot2::aes(
      step.size, aum, color=search),
      data=aum.df)+
    ggplot2::geom_line(ggplot2::aes(
      step.size, aum, color=search),
      size=1,
      data=aum.df)+
    ggplot2::facet_grid(panel ~ ., scales="free")+
    ggplot2::geom_abline(ggplot2::aes(
      slope=slope, intercept=intercept, color=search),
      data=abline.df)+
    ggplot2::geom_point(ggplot2::aes(
      0, intercept, color=search),
      data=abline.df)+
    ggplot2::scale_y_continuous("")+
    ggplot2::geom_point(
      ggplot2::aes(
        step.size, aum, color=search),
      shape=1,
      data=grid.df)
}
