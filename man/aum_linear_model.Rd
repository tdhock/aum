\name{aum_linear_model}
\alias{aum_linear_model}
\title{aum linear model}
\description{Learn a linear model with weights that minimize AUM. Weights are
initialized as a vector of zeros, then optimized using gradient
descent with exact line search.}
\usage{aum_linear_model(feature.list, 
    diff.list, max.steps = NULL, 
    improvement.thresh = NULL, 
    maxIterations = nrow(feature.list$subtrain), 
    initial.weight.fun = NULL)}
\arguments{
  \item{feature.list}{List with named elements subtrain and optionally validation, each
should be a scaled feature matrix.}
  \item{diff.list}{List with named elements subtrain and optionally validation, each
should be a data table of differences in error functions.}
  \item{max.steps}{positive integer: max number of steps of gradient descent with
exact line search (specify either this or \code{improvement.thresh}, not
both).}
  \item{improvement.thresh}{non-negative real number: keep doing gradient descent while the
improvement in AUM is greater than this number (specify either
this or \code{max.steps}, not both).}
  \item{maxIterations}{max number of iterations of exact line search, default is number
of subtrain examples.}
  \item{initial.weight.fun}{Function for computing initial weights, default NULL means use a
random standard normal vector.}
}

\value{Linear model represented as a list of class aum_linear_model with
named elements: loss is a data table of values for subtrain and
optionally validation at each step, weight.vec is the final vector
of weights learned via gradient descent, and intercept is the
value which results in minimal total error (FP+FN), learned via a
linear scan over all possible values given the final weight
vector.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}





