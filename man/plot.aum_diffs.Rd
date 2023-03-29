\name{plot.aum_diffs}
\alias{plot.aum_diffs}
\title{plot aum diffs}
\description{Plot method for \code{\link{aum_diffs}} which shows piecewise constant error
functions. Uses \code{\link{aum_errors}} internally to compute error functions
which are plotted. Not recommended for large number of examples
(>20).}
\usage{\method{plot}{aum_diffs}(x, ...)}
\arguments{
  \item{x}{data table with class "aum_diffs".}
  \item{\dots}{ignored.}
}

\value{ggplot of error functions, each example in a different panel.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}





