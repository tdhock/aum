\name{aum}
\alias{aum}
\title{aum}
\description{Compute the Area Under Minimum of False Positives and False
Negatives, and its directional derivatives.}
\usage{aum(error.diff.df, pred.vec)}
\arguments{
  \item{error.diff.df}{data frame of error differences, typically computed via
\code{\link{aum_diffs_binary}} or \code{\link{aum_diffs_penalty}}. There should be one row for
each change in error functions. "example" column indicates example
ID (int from 1 to N), "pred" column indicates predicted value
where there is a change in the error function(s), "fp_diff" and
"fn_diff" columns indicate differences in false positives and
false negatives at that predicted value. Note that this
representation assumes that each error function has fp=0 at
pred=-Inf and fn=0 at pred=Inf.}
  \item{pred.vec}{numeric vector of N predicted values.}
}

\value{Named list of two items: aum is numeric scalar loss value,
derivative_mat is N x 2 matrix of directional derivatives (first
column is derivative from left, second column is derivative from
right). If}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

(bin.diffs <- aum::aum_diffs_binary(c(0,1)))
aum::aum(bin.diffs, c(-10,10))
aum::aum(bin.diffs, c(0,0))
aum::aum(bin.diffs, c(10,-10))

}
