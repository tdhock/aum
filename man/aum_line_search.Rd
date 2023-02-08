\name{aum_line_search}
\alias{aum_line_search}
\title{aum line search}
\description{Exact line search.}
\usage{aum_line_search(error.diff.df, 
    feature.mat, weight.vec, 
    pred.vec = NULL, 
    maxIterations = nrow(error.diff.df))}
\arguments{
  \item{error.diff.df}{\code{\link{aum_diffs}} data frame with B rows, one for each breakpoint in
example-specific error functions.}
  \item{feature.mat}{N x p matrix of numeric features.}
  \item{weight.vec}{p-vector of numeric linear model coefficients.}
  \item{pred.vec}{N-vector of numeric predicted values. If NULL, \code{feature.mat} and
\code{weight.vec} will be used to compute predicted values.}
  \item{maxIterations}{positive int: max number of line search iterations.}
}

\value{List of class aum_line_search.}

\author{Toby Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

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
##         aum step.size aucAtStep aucAfterStep
## 1 0.5180666 0.0000000       0.0          0.0
## 2 0.2516401 0.1332133       0.5          1.0
## 3 0.2516401 0.4038271       0.5          0.5
## 4 0.0000000 0.5296472       1.0          1.5
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
    weight.vec=weight.vec)
  if(requireNamespace("ggplot2"))plot(nb.weight.search)

  ## Example 4: many changepoint examples, optimize predictions.
  all.ids <- rownames(neuroblastomaProcessed$feature.mat)
  all.diffs <- aum::aum_diffs_penalty(nb.err, all.ids)
  current.pred <- rep(0, length(all.ids))
  nb.all.search <- aum::aum_line_search(
    all.diffs, pred.vec=current.pred, maxIterations=2e5)
  library(data.table)
  all.result <- data.table(nb.all.search$line_search_result)
  some.result <- all.result[as.integer(seq(1, .N, l=100))]
  plot(log10(aum) ~ step.size, some.result)
  
}
