\name{aum_linear_model_cv}
\alias{aum_linear_model_cv}
\title{aum linear model cv}
\description{Cross-validation for learning number of early stopping gradient
descent steps with exact line search, in linear model for
minimizing AUM.}
\usage{aum_linear_model_cv(feature.mat, 
    diff.dt, maxIterations = "min.aum", 
    improvement.thresh = NULL, 
    n.folds = 3, initial.weight.fun = NULL)}
\arguments{
  \item{feature.mat}{N x P matrix of features, which will be scaled before gradient descent.}
  \item{diff.dt}{data table of differences in error functions, from
\code{\link{aum_diffs_penalty}} or \code{\link{aum_diffs_binary}}. There should be an example
column with values from 0 to N-1.}
  \item{maxIterations}{max iterations of the exact line search, default is number of examples.}
  \item{improvement.thresh}{before doing cross-validation to learn the number of gradient
descent steps, we do gradient descent on the full data set in
order to determine a max number of steps, by continuing to do
exact line search steps while the decrease in AUM is greater than
this value (positive real number). Default NULL means to use the
value which is ten times smaller than the min non-zero absolute
value of FP and FN diffs in \code{diff.dt}.}
  \item{n.folds}{Number of cross-validation folds to average over to determine the
best number of steps of gradient descent.}
  \item{initial.weight.fun}{Function for computing initial weight vector in gradient descent.}
}

\value{Model trained with best number of iterations, represented as a
list of class aum_linear_model_cv with named elements: keep is a
logical vector telling which features should be kept before doing
matrix multiply of learned weight vector, weight.orig/weight.vec
and intercept.orig/intercept are the learned weights/intercepts
for the original/scaled feature space, fold.loss/set.loss are data
tables of loss values for the subtrain/validation sets, used for
selecting the best number of gradient descent steps.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

## simulated binary classification problem.
N.rows <- 60
N.cols <- 2
set.seed(1)
feature.mat <- matrix(rnorm(N.rows*N.cols), N.rows, N.cols)
unknown.score <- feature.mat[,1]*2.1 + rnorm(N.rows)
label.vec <- ifelse(unknown.score > 0, 1, 0)
diffs.dt <- aum::aum_diffs_binary(label.vec)

## Default line search keeps doing iterations until increase in AUM.
(default.time <- system.time({
  default.model <- aum::aum_linear_model_cv(feature.mat, diffs.dt)
}))
plot(default.model)
print(default.valid <- default.model[["set.loss"]][set=="validation"])
print(default.model[["search"]][, .(step.size, aum, iterations=q.size)])

## Can specify max number of iterations of line search.
(small.step.time <- system.time({
  small.step.model <- aum::aum_linear_model_cv(feature.mat, diffs.dt, maxIterations = N.rows)
}))
plot(small.step.model)
print(small.step.valid <- small.step.model[["set.loss"]][set=="validation"])
small.step.model[["search"]][, .(step.size, aum, iterations=q.size)]

## Compare number of steps, iterations and time. On my machine small
## step model takes more time/steps, but less iterations in the C++
## line search code.
cbind(
  iterations=c(
    default=default.model[["search"]][, sum(q.size)],
    small.step=small.step.model[["search"]][, sum(q.size)]),
  seconds=c(
    default.time[["elapsed"]],
    small.step.time[["elapsed"]]),
  steps=c(
    default.model[["min.valid.aum"]][["step.number"]],
    small.step.model[["min.valid.aum"]][["step.number"]]),
  min.valid.aum=c(
    default.model[["min.valid.aum"]][["aum_mean"]],
    small.step.model[["min.valid.aum"]][["aum_mean"]]))

}
