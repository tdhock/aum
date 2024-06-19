\name{aum_line_search}
\alias{aum_line_search}
\title{aum line search}
\description{Exact line search using a C++ STL map (red-black tree) to
implement a queue of line intersection events. If number of rows
of \code{error.diff.df} is B, and number of iterations is I, then space
complexity is O(B) and time complexity is O( (I+B)log B ).}
\usage{aum_line_search(error.diff.df, 
    feature.mat, weight.vec, 
    pred.vec = NULL, 
    maxIterations = nrow(error.diff.df), 
    feature.mat.search = feature.mat, 
    error.diff.search = error.diff.df, 
    maxStepSize = -1)}
\arguments{
  \item{error.diff.df}{\code{\link{aum_diffs}} data frame with B rows, one for each breakpoint in
example-specific error functions.}
  \item{feature.mat}{N x p matrix of numeric features.}
  \item{weight.vec}{p-vector of numeric linear model coefficients.}
  \item{pred.vec}{N-vector of numeric predicted values. If NULL, \code{feature.mat} and
\code{weight.vec} will be used to compute predicted values.}
  \item{maxIterations}{max number of line search iterations, either a positive integer or
"max.auc" or "min.aum" indicating to keep going until AUC
decreases or AUM increases.}
  \item{feature.mat.search}{feature matrix to use in line search, default is subtrain, can be validation}
  \item{error.diff.search}{\code{\link{aum_diffs}} data frame to use in line search, default is subtrain, can be validation}
  \item{maxStepSize}{max step size to explore.}
}

\value{List of class aum_line_search. Element named "line_search_result"
is a data table with number of rows equal to \code{maxIterations} (if it
is positive integer, info for all steps, q.size column is number
of items in queue at each iteration), otherwise 1 (info for the
best step, q.size column is the total number of items popped off
the queue).}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

if(require("data.table"))setDTthreads(1L)#for CRAN check.

## Example 1: two binary data.
(bin.diffs <- aum::aum_diffs_binary(c(0,1)))
if(requireNamespace("ggplot2"))plot(bin.diffs)
bin.line.search <- aum::aum_line_search(bin.diffs, pred.vec=c(10,-10))
if(requireNamespace("ggplot2"))plot(bin.line.search)

if(requireNamespace("penaltyLearning")){

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
  X.subtrain <- X.sc[1:50,keep]
  weight.vec <- rep(0, ncol(X.subtrain))
  (diffs.subtrain <- aum::aum_diffs_penalty(nb.err, rownames(X.subtrain)))
  nb.weight.search <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec, 
    maxIterations = 200)
  if(requireNamespace("ggplot2"))plot(nb.weight.search)

  ## Stop line search after finding a (local) max AUC or min AUM.
  max.auc.search <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec,
    maxIterations="max.auc")
  min.aum.search <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec,
    maxIterations="min.aum")
  if(require("ggplot2")){
    plot(nb.weight.search)+
      geom_point(aes(
        step.size, auc),
        data=data.table(max.auc.search[["line_search_result"]], panel="auc"),
        color="red")+
      geom_point(aes(
        step.size, aum),
        data=data.table(min.aum.search[["line_search_result"]], panel="aum"),
        color="red")
  }

  ## Alternate viz with x=iteration instead of step size.
  nb.weight.full <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
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

  ## Example 4: line search on validation set.
  X.validation <- X.sc[101:300,keep]
  diffs.validation <- aum::aum_diffs_penalty(nb.err, rownames(X.validation))
  valid.search <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec, 
    maxIterations = 2000,
    feature.mat.search=X.validation,
    error.diff.search=diffs.validation)
  if(requireNamespace("ggplot2"))plot(valid.search)

  ## validation set max auc, min aum.
  max.auc.valid <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec,
    maxIterations="max.auc",
    feature.mat.search=X.validation,
    error.diff.search=diffs.validation)
  min.aum.valid <- aum::aum_line_search(
    diffs.subtrain,
    feature.mat=X.subtrain,
    weight.vec=weight.vec,
    maxIterations="min.aum",
    feature.mat.search=X.validation,
    error.diff.search=diffs.validation)
  if(require("ggplot2")){
    plot(valid.search)+
      geom_point(aes(
        step.size, auc),
        data=data.table(max.auc.valid[["line_search_result"]], panel="auc"),
        color="red")+
      geom_point(aes(
        step.size, aum),
        data=data.table(min.aum.valid[["line_search_result"]], panel="aum"),
        color="red")
  }

  ## compare subtrain and validation
  both.results <- rbind(
    data.table(valid.search$line_search_result, set="validation"),
    data.table(nb.weight.search$line_search_result, set="subtrain"))
  both.max <- rbind(
    data.table(max.auc.valid$line_search_result, set="validation"),
    data.table(max.auc.search$line_search_result, set="subtrain"))
  ggplot()+
    geom_vline(aes(
      xintercept=step.size, color=set),
      data=both.max)+
    geom_point(aes(
      step.size, auc, color=set),
      shape=1,
      data=both.results)

}

}
