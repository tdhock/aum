\name{aum_line_search_grid}
\alias{aum_line_search_grid}
\title{aum line search grid}
\description{Line search for predicted values, with grid search to check.}
\usage{aum_line_search_grid(error.diff.df, 
    feature.mat, weight.vec, 
    pred.vec = NULL, 
    maxIterations = nrow(error.diff.df), 
    n.grid = 10L, add.breakpoints = FALSE)}
\arguments{
  \item{error.diff.df}{\code{\link{aum_diffs}} data frame with B rows, one for each breakpoint in
example-specific error functions.}
  \item{feature.mat}{N x p matrix of numeric features.}
  \item{weight.vec}{p-vector of numeric linear model coefficients.}
  \item{pred.vec}{N-vector of numeric predicted values. If missing, \code{feature.mat} and
\code{weight.vec} will be used to compute predicted values.}
  \item{maxIterations}{positive int: max number of line search iterations.}
  \item{n.grid}{positive int: number of grid points for checking.}
  \item{add.breakpoints}{add breakpoints from exact search to grid search.}
}

\value{List of class aum_line_search_grid.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

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
        "\%d,\%d", intersections, intervals)),
      vjust=1.1,
      data=data.frame(
        panel="threshold", bin.line.search$line_search_result))
}

}
