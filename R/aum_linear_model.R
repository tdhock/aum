aum_linear_model_cv <- structure(function
### Cross-validation for learning number of early stopping gradient
### descent steps with exact line search, in linear model for
### minimizing AUM.
(feature.mat,
### N x P matrix of features, which will be scaled before gradient descent.
  diff.dt,
### data table of differences in error functions, from
### aum_diffs_penalty or aum_diffs_binary. There should be an example
### column with values from 0 to N-1.
  maxIterations=nrow(feature.mat),
### max iterations of the exact line search, default is number of examples.
  improvement.thresh=NULL,
### before doing cross-validation to learn the number of gradient
### descent steps, we do gradient descent on the full data set in
### order to determine a max number of steps, by continuing to do
### exact line search steps while the decrease in AUM is greater than
### this value (positive real number). Default NULL means to use the
### value which is ten times smaller than the min non-zero absolute
### value of FP and FN diffs in diff.dt.
  n.folds=3,
### Number of cross-validation folds to average over to determine the
### best number of steps of gradient descent.
  initial.weight.fun=NULL
### Function for computing initial weight vector in gradient descent.
){
  . <- fp_diff <- fn_diff <- example <- fp <- fn <- fold <- pred <- 
    valid.fold <- sd <- aum_mean <- NULL
  ## Above to suppress CRAN NOTE.
  example.totals <- diff.dt[, .(
    fn=sum(fn_diff),
    fp=sum(fp_diff)
  ), by=example]
  
  if(is.null(improvement.thresh)){
    abs.diff <- diff.dt[, abs(c(fp_diff, fn_diff))]
    not.zero <- abs.diff[0 < abs.diff]
    improvement.thresh <- min(not.zero)/10
    ## TODO: does this heuristic generalize well to other data sets?
  }
  X.sc <- scale(feature.mat)
  keep <- apply(is.finite(X.sc), 2, all)
  X.keep <- X.sc[,keep,drop=FALSE]
  train.features <- list(subtrain=X.keep)
  train.diffs <- list(subtrain=diff.dt)
  overfit.model <- aum_linear_model(
    train.features, train.diffs,
    initial.weight.fun=initial.weight.fun,
    improvement.thresh=improvement.thresh,
    maxIterations=maxIterations)
  uniq.folds <- 1:n.folds
  zero.counts <- colSums(example.totals[, .(fn,fp)]==0)
  minority <- names(zero.counts)[which.max(zero.counts)]
  minority.zero <- example.totals[[minority]]==0
  example.totals[, fold := sample(
    rep(sample(uniq.folds), l=.N)
  ), by=minority.zero]
  minority.folds <- example.totals[minority.zero==FALSE, length(unique(fold))]
  if(minority.folds < n.folds){
    stop(sprintf("not enough data for %d-fold cross-validation, because there are only %d examples for which there are non-zero values for the minority diff, %s", n.folds, minority.folds, minority))
  }
  fold.loss <- data.table(valid.fold=uniq.folds)[, {
    logical.list <- with(example.totals, list(
      subtrain=fold!=valid.fold,
      validation=fold==valid.fold))
    diff.list <- lapply(logical.list, function(is.set){
      some.indices <- which(is.set)
      all.indices <- rep(NA, nrow(X.keep))
      all.indices[some.indices] <- seq_along(some.indices)-1L
      diff.dt[, .(
        example=all.indices[example+1L], pred, fp_diff, fn_diff
      )][!is.na(example)]
    })
    feature.list <- lapply(logical.list, function(is.set){
      X.keep[is.set,]
    })
    valid.model <- aum_linear_model(
      feature.list, diff.list,
      initial.weight.fun=initial.weight.fun,
      max.steps=max(overfit.model$loss$step.number),
      maxIterations=maxIterations)
    valid.model$loss
  }, by=valid.fold]
  set.loss <- dcast(
    fold.loss,
    step.number + set ~ .,
    list(mean, sd),
    value.var="aum")
  best.row <- set.loss[set=="validation"][which.min(aum_mean)]
  final.model <- aum_linear_model(
    train.features, train.diffs,
    initial.weight.fun=initial.weight.fun,
    max.steps=best.row$step.number,
    maxIterations=maxIterations)
  final.model$fold.loss <- fold.loss
  final.model$set.loss <- set.loss
  final.model$keep <- keep
  final.model$weight.orig <-
    final.model$weight.vec/attr(X.sc, "scaled:scale")[keep]
  final.model$intercept.orig <- final.model$intercept-sum(
    final.model$weight.orig*attr(X.sc, "scaled:center")[keep])
  structure(final.model, class="aum_linear_model_cv")
### Model trained with best number of iterations, represented as a
### list of class aum_linear_model_cv with named elements: keep is a
### logical vector telling which features should be kept before doing
### matrix multiply of learned weight vector, weight.orig/weight.vec
### and intercept.orig/intercept are the learned weights/intercepts
### for the original/scaled feature space, fold.loss/set.loss are data
### tables of loss values for the subtrain/validation sets, used for
### selecting the best number of gradient descent steps.
}, ex=function(){

  ## simulated binary classification problem.
  N.rows <- 100
  N.cols <- 20
  set.seed(1)
  feature.mat <- matrix(rnorm(N.rows*N.cols), N.rows, N.cols)
  unknown.score <- feature.mat[,1]*2.1 + rnorm(N.rows)
  label.vec <- ifelse(unknown.score > 3, 1, 0)
  diffs.dt <- aum::aum_diffs_binary(label.vec)
  model <- aum::aum_linear_model_cv(feature.mat, diffs.dt)
  plot(model)
  
})

predict.aum_linear_model_cv <- function(object, newdata, ...){
  newdata[,object$keep,drop=FALSE] %*% object$weight.orig +
    object$intercept.orig
}

plot.aum_linear_model_cv <- function(x, ...){
  lattice::xyplot(
    aum_mean ~ step.number, x$set.loss, 
    groups=set, type="l", 
    auto.key=list(space="right", points=FALSE, lines=TRUE))
}

aum_linear_model <- function
### Learn a linear model with weights that minimize AUM. Weights are
### initialized as a vector of zeros, then optimized using gradient
### descent with exact line search.
(feature.list,
### List with named elements subtrain and optionally validation, each
### should be a scaled feature matrix.
  diff.list,
### List with named elements subtrain and optionally validation, each
### should be a data table of differences in error functions.
  max.steps=NULL,
### positive integer: max number of steps of gradient descent with
### exact line search (specify either this or improvement.thresh, not
### both).
  improvement.thresh=NULL,
### non-negative real number: keep doing gradient descent while the
### improvement in AUM is greater than this number (specify either
### this or max.steps, not both).
  maxIterations=nrow(feature.list$subtrain),
### max number of iterations of exact line search, default is number
### of subtrain examples.
  initial.weight.fun=NULL
### Function for computing initial weights, default NULL means use a
### random standard normal vector.
){
  fp_before <- fn_before <- thresh <- NULL
  ## Above to suppress CRAN NOTE.
  weight.vec <- if(is.null(initial.weight.fun)){
    rnorm(ncol(feature.list$subtrain))
  }else{
    initial.weight.fun(feature.list$subtrain, diff.list$subtrain)
  }
  improvement <- old.aum <- Inf
  step.number <- 0
  loss.dt.list <- list()
  while({
    search.result <- aum::aum_line_search(
      diff.list$subtrain,
      maxIterations=maxIterations,
      feature.mat=feature.list$subtrain,
      weight.vec=weight.vec)
    loss.dt.list[[paste(step.number, "subtrain")]] <- data.table(
      step.number, 
      set="subtrain",
      aum=search.result$aum)
    if("validation"%in%names(feature.list)){
      valid.list <- aum::aum(
        diff.list$validation,
        feature.list$validation %*% weight.vec)
      loss.dt.list[[paste(step.number, "validation")]] <- data.table(
        step.number,
        set="validation",
        aum=valid.list$aum)
    }
    exact.dt <- data.table(search.result$line_search_result)
    best.row <- exact.dt[which.min(aum)]
    improvement <- old.aum-best.row$aum
    old.aum <- best.row$aum
    if(!is.null(improvement.thresh)){
      improvement.thresh < improvement
    }else if(!is.null(max.steps)){
      step.number < max.steps 
    }else{
      stop("either improvement.thresh or max.steps must be not NULL")
    }
  }){
    step.number <- step.number+1
    weight.vec <- weight.vec-
      best.row$step.size*search.result$gradient_weight
  }
  out.list <- list(
    loss=do.call(rbind, loss.dt.list),
    weight.vec=weight.vec,
    intercept=data.table(
      search.result$total_error, key="thresh"
    )[,{
      best <- which.min(fp_before+fn_before)
      if(best==1){
        thresh[1]-1
      }else{
        mean(thresh[c(best-1,best)])
      }
    }])
  structure(out.list, class="aum_linear_model")
### Linear model represented as a list of class aum_linear_model with
### named elements: loss is a data table of values for subtrain and
### optionally validation at each step, weight.vec is the final vector
### of weights learned via gradient descent, and intercept is the
### value which results in minimal total error (FP+FN), learned via a
### linear scan over all possible values given the final weight
### vector.
}  

plot.aum_linear_model <- function(x, ...){
  lattice::xyplot(
    aum ~ step.number, x$loss, 
    groups=set, type="l", 
    auto.key=list(space="right", points=FALSE, lines=TRUE))
}
