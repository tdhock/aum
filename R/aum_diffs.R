aum_errors <- structure(function
### Convert diffs to canonical errors, used internally in
### plot.aum_diffs.
(diffs.df
### data.table of diffs from aum_diffs.
){
  pred <- fp_diff <- fn_diff <- example <- . <- NULL
  ## Above to silence CRAN check NOTE.
  data.table(diffs.df)[order(example, pred), .(
    min.pred=c(-Inf, pred),
    max.pred=c(pred, Inf),
    fp=cumsum(c(0, fp_diff)),
    fn=rev(cumsum(c(0, -rev(fn_diff))))
  ), by=example]
### data.table suitable for plotting piecewise constant error
### functions, with columns example, min.pred, max.pred, fp, fn.
}, ex=function(){

  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  if(requireNamespace("ggplot2"))plot(bin.diffs)
  aum::aum_errors(bin.diffs)

})

plot.aum_diffs <- function
### Plot method for aum_diffs which shows piecewise constant error
### functions. Uses aum_errors internally to compute error functions
### which are plotted. Not recommended for large number of examples
### (>20).
(x,
### data table with class "aum_diffs".
  ...
### ignored.
){
  min.pred <- value <- max.pred <- variable <- pred <- NULL
  ## Above to silence CRAN check NOTE.
  err.wide <- aum_errors(x)
  err.tall <- data.table::melt(err.wide, measure=c("fp", "fn"))
  ggplot2::ggplot()+
    ggplot2::geom_segment(ggplot2::aes(
      min.pred, value,
      xend=max.pred, yend=value,
      color=variable, size=variable),
      data=err.tall)+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=pred),
      data=x)+
    ggplot2::scale_size_manual(values=c(fp=2, fn=1))+
    ggplot2::facet_grid(example ~ ., labeller=ggplot2::label_both)
### ggplot of error functions, each example in a different panel.
}

aum_diffs <- structure(function
### Create error differences data table which can be used as input to
### aum function. Typical users should not use this function directly,
### and instead use aum_diffs_binary for binary classification, and
### aum_diffs_penalty for error defined as a function of non-negative
### penalty.
(example,
### Integer or character vector identifying different examples.
  pred,
### Numeric vector of predicted values at which the error changes.
  fp_diff,
### Numeric vector of difference in fp at pred.
  fn_diff,
### Numeric vector of difference in fn at pred.
  pred.name.vec
### Character vector of example names for predictions.
){
  if(is.character(example) && is.character(pred.name.vec)){
    n.tab <- table(pred.name.vec)
    bad <- n.tab[1 < n.tab]
    if(length(bad)){
      stop(
        "elements of pred.name.vec must be unique, problems: ",
        paste(names(bad), collapse=", "))
    }
    name2id <- seq(0, length(pred.name.vec)-1L)
    names(name2id) <- pred.name.vec
    example <- name2id[example]
  }
  if(!is.integer(example)){
    stop("example must be integer vector but has class: ", class(example)[1])
  }
  out <- data.table(example, pred, fp_diff, fn_diff)[!is.na(example)]
  setkey(out, example, pred)
  class(out) <- c("aum_diffs", class(out))
  out
### data table of class "aum_diffs" in which each rows represents a
### breakpoint in an error function. Columns are interpreted as
### follows: there is a change of "fp_diff","fn_diff" at predicted
### value "pred" for example/observation "example". This can be used
### for computing Area Under Minimum via aum function, and plotted via
### plot.aum_diffs.
}, ex=function(){

  aum::aum_diffs_binary(c(0,1))
  aum::aum_diffs(c("positive", "negative"), 0, c(0,1), c(-1,1), c("negative", "positive"))
  rbind(aum::aum_diffs(0L, 0, 1, 0), aum_diffs(1L, 0, 0, -1))

})

aum_diffs_binary <- structure(function
### Convert binary labels to error differences.
(label.vec,
### Numeric vector representing binary labels (either all 0,1 or all
### -1,1). If named, names are used to identify each example.
  pred.name.vec,
### Character vector of prediction example names, used to convert
### names of label.vec to integers.
  denominator="count"
### Type of diffs, either "count" or "rate".
){
  allin <- function(...)all(label.vec %in% c(...))
  if(!all(
    is.numeric(label.vec),
    0<length(label.vec),
    is.finite(label.vec),
    allin(0,1) || allin(-1,1)
  )){
    stop("label.vec must be numeric vector with length>0 and all elements either 0,1 or -1,1")
  }
  is.positive <- label.vec==1
  example <- if(is.null(names(label.vec))){
    seq(0, length(label.vec)-1)
  }else{
    n.vec <- names(label.vec)
    n.tab <- table(n.vec)
    bad <- n.tab[1 < n.tab]
    if(length(bad)){
      stop(
        "if label.vec has names they must be unique, problems: ",
        paste(names(bad), collapse=", "))
    }
    n.vec
  }
  if(identical(denominator, "rate")){
    fp.denom <- sum(!is.positive)
    fn.denom <- sum(is.positive)
  }else{
    fp.denom <- 1
    fn.denom <- 1
  }
  aum_diffs(
    example,
    0,
    ifelse(is.positive,  0, 1)/fp.denom,
    ifelse(is.positive, -1, 0)/fn.denom,
    pred.name.vec)
### data table of class "aum_diffs" in which each rows represents a
### breakpoint in an error function. Columns are interpreted as
### follows: there is a change of "fp_diff","fn_diff" at predicted
### value "pred" for example/observation "example". This can be used
### for computing Area Under Minimum via aum function, and plotted via
### plot.aum_diffs.
}, ex=function(){

  aum_diffs_binary(c(0,1))
  aum_diffs_binary(c(-1,1))
  aum_diffs_binary(c(a=0,b=1,c=0), pred.name.vec=c("c","b"))
  aum_diffs_binary(c(0,0,1,1,1), denominator="rate")

})

aum_diffs_penalty <- structure(function
### Convert penalized errors to error differences. A typical use case
### is for penalized optimal changepoint models, for which small
### penalty values result in large fp/fn, and large penalty values
### result in small fp/fn.
(errors.df,
### data.frame which describes error as a function of penalty/lambda,
### with at least columns example, min.lambda, fp, fn. Interpreted as
### follows: fp/fn occur from all penalties from min.lambda to the
### next value of min.lambda within the current value of example.
  pred.name.vec,
### Character vector of prediction example names, used to convert
### names of label.vec to integers.
  denominator="count"
### Type of diffs, either "count" or "rate".
){
  example <- min.lambda <- fp <- fn <- n.zero <- more <- . <- NULL
  ## Above to silence CRAN check NOTE.
  for(cname in c("fp", "fn", "min.lambda")){
    if(!is.numeric(errors.df[[cname]])){
      stop("errors.df must have numeric column named ", cname)
    }
  }
  e <- errors.df[["example"]]
  if(!(is.integer(e) || is.character(e))){
    stop("errors.df must have integer or character column named example")
  }
  err.dt <- as.data.table(errors.df)[order(example, -min.lambda)]
  if(!missing(pred.name.vec)){
    err.dt <- err.dt[example %in% pred.name.vec]
  }
  prob.dt <- err.dt[, {
    tab <- table(min.lambda)
    data.table(
      n.zero=sum(min.lambda==0),
      more=paste(names(tab)[1 < tab], collapse=","))
  }, by=example]
  no.zero <- prob.dt[n.zero < 1]
  if(nrow(no.zero)){
    stop(
      "need min.lambda=0 for each example, problems: ",
      paste(no.zero$example, collapse=","))
  }
  more.dt <- prob.dt[more != ""]
  if(nrow(more.dt)){
    bad.ex.vec <- more.dt[, paste0(example, ":", more)]
    stop(
      "need only one min.lambda per example, problems with more are (example:min.lambda) ",
      paste(bad.ex.vec, collapse=" "))
  }
  with(err.dt, {
    is.end <- min.lambda == 0
    mydiff <- function(x){
      ifelse(is.end, 0, diff(x))
    }
    fp_diff <- mydiff(fp)
    fn_diff <- mydiff(fn)
    keep <- fp_diff != 0 | fn_diff != 0
    if(identical(denominator, "rate")){
      fp.denom <- sum(fp_diff)
      fn.denom <- -sum(fn_diff)
    }else{
      fp.denom <- 1
      fn.denom <- 1
    }
    aum_diffs(
      example[keep],
      -log(min.lambda[keep]),
      fp_diff[keep]/fp.denom,
      fn_diff[keep]/fn.denom,
      pred.name.vec)
  })
### data table of class "aum_diffs" in which each rows represents a
### breakpoint in an error function. Columns are interpreted as
### follows: there is a change of "fp_diff","fn_diff" at predicted
### value "pred" for example/observation "example". This can be used
### for computing Area Under Minimum via aum function, and plotted via
### plot.aum_diffs.
}, ex=function(){

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

})

