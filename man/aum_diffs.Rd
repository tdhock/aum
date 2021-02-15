\name{aum_diffs}
\alias{aum_diffs}
\title{aum diffs}
\description{Create error differences data table which can be used as input to
\code{\link{aum}} function.}
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

\value{data table with "aum_diffs" class and same columns as input
arguments.}

\author{Toby Dylan Hocking}




\examples{

aum::aum_diffs_binary(c(0,1))
aum::aum_diffs(c("positive", "negative"), 0, c(0,1), c(-1,1), c("negative", "positive"))
rbind(aum::aum_diffs(0L, 0, 1, 0), aum_diffs(1L, 0, 0, -1))

}
