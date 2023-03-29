\name{aum_diffs}
\alias{aum_diffs}
\title{aum diffs}
\description{Create error differences data table which can be used as input to
\code{\link{aum}} function. Typical users should not use this function directly,
and instead use \code{\link{aum_diffs_binary}} for binary classification, and
\code{\link{aum_diffs_penalty}} for error defined as a function of non-negative
penalty.}
\usage{aum_diffs(example, pred, 
    fp_diff, fn_diff, 
    pred.name.vec)}
\arguments{
  \item{example}{Integer or character vector identifying different examples.}
  \item{pred}{Numeric vector of predicted values at which the error changes.}
  \item{fp_diff}{Numeric vector of difference in fp at \code{pred}.}
  \item{fn_diff}{Numeric vector of difference in fn at \code{pred}.}
  \item{pred.name.vec}{Character vector of \code{example} names for predictions.}
}

\value{data table of class "aum_diffs" in which each rows represents a
breakpoint in an error function. Columns are interpreted as
follows: there is a change of "fp_diff","fn_diff" at predicted
value "pred" for example/observation "example". This can be used
for computing Area Under Minimum via \code{\link{aum}} function, and plotted via
\code{\link{plot.aum_diffs}}.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

aum::aum_diffs_binary(c(0,1))
aum::aum_diffs(c("positive", "negative"), 0, c(0,1), c(-1,1), c("negative", "positive"))
rbind(aum::aum_diffs(0L, 0, 1, 0), aum_diffs(1L, 0, 0, -1))

}
