aum_errors <- structure(function
### Convert diffs to canonical errors.
(diffs.dt
### data.table of diffs from aum_diffs.
){
  pred <- fp_diff <- fn_diff <- example <- NULL
  ## Above to silence CRAN check NOTE.
  diffs.dt[, data.table(
    min.pred=c(-Inf, pred),
    max.pred=c(pred, Inf),
    fp=cumsum(c(0, fp_diff)),
    fn=rev(cumsum(c(0, -rev(fn_diff))))
  ), by=example]
### data.table suitable for plotting piecewise constant error
### functions, with columns example, min.pred, max.pred, fp, fn.
}, ex=function(){

  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  plot(bin.diffs)
  aum::aum_errors(bin.diffs)
  
})

### Plot method for aum_diffs which shows piecewise constant error
### functions.
plot.aum_diffs <- function(x, ...){
  min.pred <- value <- max.pred <- variable <- pred <- NULL
  ## Above to silence CRAN check NOTE.
  if(!requireNamespace("ggplot2")){
    stop("please install ggplot2 for plotting aum_diffs")
  }
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
}  

aum_diffs <- structure(function
### Create error differences data table which can be used as input to
### aum function.
(example,
### Integer or character vector identifying different examples.
  pred,
### Numeric vector of predicted values at which the error changes.
  fp_diff,
### Numeric vector of difference in fp at pred.
  fn_diff
### Numeric vector of difference in fn at pred.
){
  out <- data.table(example, pred, fp_diff, fn_diff)
  class(out) <- c("aum_diffs", class(out))
  out
### data table with "aum_diffs" class and same columns as input
### arguments.
}, ex=function(){

  aum::aum_diffs_binary(c(0,1))
  rbind(aum::aum_diffs(1, 0, 1, 0), aum_diffs(2, 0, 0, -1))
  
})

aum_diffs_binary <- structure(function
### Convert binary labels to error differences.
(label.vec
### Numeric vector representing binary labels (either all 0,1 or all
### -1,1).
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
  label.vec[label.vec==0] <- -1
  example <- if(is.null(names(label.vec))){
    seq_along(label.vec)
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
  aum_diffs(
    example,
    0,
    ifelse(label.vec==1,  0, 1),
    ifelse(label.vec==1, -1, 0))
### data.frame of error diffs which can be used as input to the aum
### function.
}, ex=function(){

  aum_diffs_binary(c(0,1))
  aum_diffs_binary(c(-1,1))
  aum_diffs_binary(c(a=0,b=1,c=0))
  
})

aum_diffs_penalty <- structure(function
### Convert penalized errors to error differences. A typical use case
### is for penalized optimal changepoint models, for which small
### penalty values result in large fp/fn, and large penalty values
### result in small fp/fn.
(errors.df
### data.frame which describes error as a function of penalty/lambda,
### with at least columns example, min.lambda, fp, fn. Interpreted as
### follows: fp/fn occur from all penalties from min.lambda to the
### next value of min.lambda within the current value of example.
){
  example <- min.lambda <- fp <- fn <- NULL
  ## Above to silence CRAN check NOTE.
  with(as.data.table(errors.df)[order(example, -min.lambda)], {
    is.end <- min.lambda == 0
    mydiff <- function(x){
      ifelse(is.end, 0, diff(x))
    }
    fp_diff <- mydiff(fp)
    fn_diff <- mydiff(fn)
    keep <- fp_diff != 0 | fn_diff != 0
    aum_diffs(
      example[keep],
      -log(min.lambda[keep]),
      fp_diff[keep],
      fn_diff[keep])
  })
### data table of error diffs which can be used as input to the aum
### function.
}, ex=function(){

  ## Simple synthetic example with two changes in error function.
  simple.df <- data.frame(
    example=1,
    min.lambda=c(0, exp(1), exp(2), exp(3)),
    fp=c(6,2,2,0),
    fn=c(0,1,1,5))
  (simple.diffs <- aum::aum_diffs_penalty(simple.df))
  plot(simple.diffs)

  ## Simple real data with four example, one has non-monotonic fn.
  if(requireNamespace("penaltyLearning")){
    data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
    ## assume min.lambda, max.lambda columns only? use names?
    nb.err <- with(neuroblastomaProcessed$errors, data.frame(
      example=paste0(profile.id, ".", chromosome),
      min.lambda,
      max.lambda,
      fp, fn))
    nb.some <- subset(nb.err, example %in% c("1.1", "1.2", "4.1", "4.2"))
    aum::aum_diffs_penalty(nb.some)
  }

  ## More complex real data example
  data(fn.not.zero, package="aum", envir=environment())
  (fn.not.zero.diffs <- aum::aum_diffs_penalty(fn.not.zero))

  if(require("ggplot2")){

    ex <- function(x)gsub("/", "\n", x)
    fn.not.zero[, ex := ex(example)]
    fn.not.zero.diffs[, ex := ex(example)]
    fn.not.zero.tall <- data.table::melt(fn.not.zero, measure=c("fp", "fn"))
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
      facet_grid(ex ~ .)

  }
  
})
  
