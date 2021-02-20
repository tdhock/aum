\name{aum_diffs_binary}
\alias{aum_diffs_binary}
\title{aum diffs binary}
\description{Convert binary labels to error differences.}
\usage{aum_diffs_binary(label.vec, 
    pred.name.vec, denominator = "count")}
\arguments{
  \item{label.vec}{Numeric vector representing binary labels (either all 0,1 or all
-1,1). If named, names are used to identify each example.}
  \item{pred.name.vec}{Character vector of prediction example names, used to convert
names of \code{label.vec} to integers.}
  \item{denominator}{Type of diffs, either "count" or "rate".}
}

\value{data.frame of error diffs which can be used as input to the \code{\link{aum}}
function.}

\author{Toby Dylan Hocking}




\examples{

aum_diffs_binary(c(0,1))
aum_diffs_binary(c(-1,1))
aum_diffs_binary(c(a=0,b=1,c=0), pred.name.vec=c("c","b"))
aum_diffs_binary(c(0,0,1,1,1), denominator="rate")

}
