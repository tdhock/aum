\name{fn.not.zero}
\alias{fn.not.zero}
\docType{data}
\title{
  Penalized models with non-zero fn at penalty=0
}
\description{
  Usually we assume that fn must be zero at penalty=0, but this is not
  always the case in real data/labels. For example in the PeakSegDisk
  model with penalty=0, there are peaks almost everywhere but if a positive
  label is too small or misplaced with respect to the detected peaks,
  then there can be false negatives.
}
\usage{data("fn.not.zero")}
\format{
  A data frame with 156 observations on the following 5 variables.
  \describe{
    \item{\code{example}}{a character vector}
    \item{\code{min.lambda}}{a numeric vector}
    \item{\code{max.lambda}}{a numeric vector}
    \item{\code{fp}}{a numeric vector}
    \item{\code{fn}}{a numeric vector}
  }
}
\source{
  https://github.com/tdhock/feature-learning-benchmark
}
\keyword{datasets}
