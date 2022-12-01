aum_line_search <- structure(function
(error.diff.df,
  pred.vec,
  maxIterations=nrow(error.diff.df)){
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
  L
}, ex=function(){

  ## Example 1: two binary data.
  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  aum::aum_line_search(bin.diffs, c(10,-10))

  ## Example 2: two changepoint examples, one with three breakpoints.
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  (nb.diffs <- aum::aum_diffs_penalty(nb.err, c("1.1", "4.2")))
  if(requireNamespace("ggplot2"))plot(nb.diffs)
  aum::aum_line_search(nb.diffs, c(10,-10))
  
})
