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

## Example 5: all changepoint examples, with noise features and linear model.
signal.features <- neuroblastomaProcessed$feature.mat[,c("log2.n","log.hall")]
n.noise <- 20
set.seed(1)
noise.features <- matrix(
  rnorm(n.noise*nrow(signal.features)),
  nrow(signal.features), n.noise)
X.sc <- scale(cbind(signal.features, noise.features))
keep <- apply(is.finite(X.sc), 2, all)
X.keep <- X.sc[,keep]
subtrain.i <- sample(1:nrow(X.keep), nrow(X.keep)/2)
index.list <- list(
  subtrain=subtrain.i,
  validation=-subtrain.i)
diff.list <- lapply(index.list, function(set.i){
  aum::aum_diffs_penalty(nb.err, rownames(X.keep)[set.i])
})
weight.vec <- rep(0, ncol(X.keep))
improvement <- old.aum <- Inf
iteration <- 0
loss.dt.list <- list()
while(improvement > 1e-4){
  iteration <- iteration+1
  valid.list <- aum::aum(
    diff.list$validation, X.keep[index.list$validation,] \%*\% weight.vec)
  nb.weight.search <- aum::aum_line_search(
    diff.list$subtrain,
    feature.mat=X.keep[index.list$subtrain,],
    weight.vec=weight.vec)
  loss.dt.list[[iteration]] <- data.table(
    iteration,
    set=c("subtrain", "validation"),
    aum=c(nb.weight.search$aum, valid.list$aum))
  exact.dt <- data.table(nb.weight.search$line_search_result)
  exact.dt[, kink := .I/.N]
  best.row <- exact.dt[which.min(aum)]
  if(interactive())cat(sprintf(
    "iteration=\%4d aum=\%.6f step=\%.6f kink=\%f\\n",
    iteration, best.row$aum, best.row$step.size, best.row$kink))
  improvement <- old.aum-best.row$aum
  old.aum <- best.row$aum
  weight.vec <- weight.vec-best.row$step.size*nb.weight.search$gradient_weight
}
loss.dt <- do.call(rbind, loss.dt.list)
min.dt <- loss.dt[, .SD[which.min(aum)], by=set]
if(requireNamespace("ggplot2")){
  ggplot2::ggplot()+
    ggplot2::geom_line(ggplot2::aes(
      iteration, aum, color=set),
      data=loss.dt)+
    ggplot2::geom_point(ggplot2::aes(
      iteration, aum, color=set),
      shape=1,
      data=min.dt)+
    ggplot2::scale_y_log10()
}

}
