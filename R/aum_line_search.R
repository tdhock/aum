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
  pred.i <- error.diff.df$example+1L
  L$line_search_input <- data.frame(
    fp.diff=error.diff.df$fp_diff,
    fn.diff=error.diff.df$fn_diff,
    intercept=error.diff.df$pred-pred.vec[pred.i],
    slope=rowMeans(L$derivative_mat)[pred.i])
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
  aum.df <- data.frame(panel="aum", x$line_search_result)
  abline.df <- data.frame(panel="threshold", x$line_search_input)
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0,"lines"))+
    geom_vline(aes(
      xintercept=step.size),
      color="grey",
      data=x$line_search_result)+
    geom_point(aes(
      step.size, aum),
      data=aum.df)+
    geom_line(aes(
      step.size, aum),
      size=1,
      data=aum.df)+
    facet_grid(panel ~ ., scales="free")+
    geom_abline(aes(
      slope=slope, intercept=intercept),
      size=1,
      data=abline.df)+
    geom_point(aes(
      0, intercept),
      data=abline.df)+
    scale_y_continuous("")
}
