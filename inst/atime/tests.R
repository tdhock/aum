data(neuroblastomaProcessed, package="penaltyLearning")
nb.err <- with(neuroblastomaProcessed$errors, data.frame(
  example=paste0(profile.id, ".", chromosome),
  min.lambda,
  max.lambda,
  fp, fn))
X.sc <- scale(neuroblastomaProcessed$feature.mat)
keep <- apply(is.finite(X.sc), 2, all)
X.keep <- X.sc[,keep]
weight.vec <- rep(0, ncol(X.keep))
nb.diffs <- aum::aum_diffs_penalty(nb.err, rownames(X.keep))

biggest.max.it <- 4278103
logseq <- function(x,by=0.25){
  as.integer(10^seq(1, log10(x), by=by))
}
test.list <- list(
  "N=data,iterations=N"=list(
    expr=quote(aum::aum_line_search(
      some.diffs,
      feature.mat=X.keep,
      weight.vec=weight.vec, 
      maxIterations = nrow(some.diffs))),
    N=logseq(nrow(nb.diffs)),
    setup=quote(some.diffs <- nb.diffs[1:N])
  ),
  "N=iterations,data=full"=list(
    expr=quote(aum::aum_line_search(
      nb.diffs,
      feature.mat=X.keep,
      weight.vec=weight.vec, 
      maxIterations = N)),
    N=logseq(biggest.max.it),
    setup={}
  )
)
