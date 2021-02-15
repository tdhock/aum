\name{aum_errors}
\alias{aum_errors}
\title{aum errors}
\description{Convert diffs to canonical errors.}
\usage{aum_errors(diffs.dt)}
\arguments{
  \item{diffs.dt}{data.table of diffs from \code{\link{aum_diffs}}.}
}

\value{data.table suitable for plotting piecewise constant error
functions, with columns example, min.pred, max.pred, fp, fn.}

\author{Toby Dylan Hocking}




\examples{

(bin.diffs <- aum::aum_diffs_binary(c(0,1)))
if(requireNamespace("ggplot2"))plot(bin.diffs)
aum::aum_errors(bin.diffs)

}
