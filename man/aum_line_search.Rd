\name{aum_line_search}
\alias{aum_line_search}
\title{aum line search}
\description{Line search for predicted values.}
\usage{aum_line_search(error.diff.df, 
    pred.vec, maxIterations = nrow(error.diff.df))}
\arguments{
  \item{error.diff.df}{\code{\link{aum_diffs}} data frame.}
  \item{pred.vec}{Numeric vector of predicted values.}
  \item{maxIterations}{positive int: max number of line search iterations.}
}

\value{List of class aum_line_search.}

\author{Toby Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

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

}
