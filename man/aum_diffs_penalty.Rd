\name{aum_diffs_penalty}
\alias{aum_diffs_penalty}
\title{aum diffs penalty}
\description{Convert penalized errors to error differences. A typical use case
is for penalized optimal changepoint models, for which small
penalty values result in large fp/fn, and large penalty values
result in small fp/fn.}
\usage{aum_diffs_penalty(errors.df, 
    pred.name.vec, denominator = "count")}
\arguments{
  \item{errors.df}{data.frame which describes error as a function of penalty/lambda,
with at least columns example, min.lambda, fp, fn. Interpreted as
follows: fp/fn occur from all penalties from min.lambda to the
next value of min.lambda within the current value of example.}
  \item{pred.name.vec}{Character vector of prediction example names, used to convert
names of label.vec to integers.}
  \item{denominator}{Type of diffs, either "count" or "rate".}
}

\value{data table of class "aum_diffs" in which each rows represents a
breakpoint in an error function. Columns are interpreted as
follows: there is a change of "fp_diff","fn_diff" at predicted
value "pred" for example/observation "example". This can be used
for computing Area Under Minimum via \code{\link{aum}} function, and plotted via
\code{\link{plot.aum_diffs}}.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

## Simple synthetic example with two changes in error function.
simple.df <- data.frame(
  example=1L,
  min.lambda=c(0, exp(1), exp(2), exp(3)),
  fp=c(6,2,2,0),
  fn=c(0,1,1,5))
(simple.diffs <- aum::aum_diffs_penalty(simple.df))
if(requireNamespace("ggplot2"))plot(simple.diffs)
(simple.rates <- aum::aum_diffs_penalty(simple.df, denominator="rate"))
if(requireNamespace("ggplot2"))plot(simple.rates)

## Simple real data with four example, one has non-monotonic fn.
if(requireNamespace("penaltyLearning")){
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  ## assume min.lambda, max.lambda columns only? use names?
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  (nb.diffs <- aum::aum_diffs_penalty(nb.err, c("1.2", "1.1", "4.1", "4.2")))
  if(requireNamespace("ggplot2"))plot(nb.diffs)
}

## More complex real data example
data(fn.not.zero, package="aum", envir=environment())
pred.names <- unique(fn.not.zero$example)
(fn.not.zero.diffs <- aum::aum_diffs_penalty(fn.not.zero, pred.names))
if(requireNamespace("ggplot2"))plot(fn.not.zero.diffs)

if(require("ggplot2")){
  name2id <- structure(seq(0, length(pred.names)-1L), names=pred.names)
  fn.not.zero.wide <- fn.not.zero[, .(example=name2id[example], min.lambda, max.lambda, fp, fn)]
  fn.not.zero.tall <- data.table::melt(fn.not.zero.wide, measure=c("fp", "fn"))
  ggplot()+
    geom_segment(aes(
      -log(min.lambda), value,
      xend=-log(max.lambda), yend=value,
      color=variable, size=variable),
      data=fn.not.zero.tall)+
    geom_point(aes(
      -log(min.lambda), value,
      fill=variable),
      color="black",
      shape=21,
      data=fn.not.zero.tall)+
    geom_vline(aes(
      xintercept=pred),
      data=fn.not.zero.diffs)+
    scale_size_manual(values=c(fp=2, fn=1))+
    facet_grid(example ~ ., labeller=label_both)
}

}
