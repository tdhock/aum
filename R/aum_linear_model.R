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
  n.folds=3
### Number of cross-validation folds to average over to determine the
### best number of steps of gradient descent.
){
  if(is.null(improvement.thresh)){
    abs.diff <- diff.dt[, abs(c(fp_diff, fn_diff))]
    not.zero <- abs.diff[0 < abs.diff]
    improvement.thresh <- min(not.zero)/10
  }
  X.sc <- scale(feature.mat)
  keep <- apply(is.finite(X.sc), 2, all)
  X.keep <- X.sc[,keep,drop=FALSE]
  train.features <- list(subtrain=X.keep)
  train.diffs <- list(subtrain=diff.dt)
  overfit.model <- aum_linear_model(
    train.features, train.diffs,
    improvement.thresh=improvement.thresh,
    maxIterations=maxIterations)
  uniq.folds <- 1:n.folds
  fold.vec <- sample(rep(uniq.folds, l=nrow(X.keep)))
  fold.loss <- data.table(valid.fold=uniq.folds)[, {
    logical.list <- list(
      subtrain=fold.vec!=valid.fold,
      validation=fold.vec==valid.fold)
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

  ## learn a model for a real changepoint data set.
  if(requireNamespace("penaltyLearning")){

    library(data.table)
    data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
    nb.err <- with(neuroblastomaProcessed$errors, data.table(
      example=paste0(profile.id, ".", chromosome),
      min.log.lambda,
      max.log.lambda, 
      min.lambda, 
      max.lambda, 
      fp, fn, errors=fp+fn))
    signal.features <- neuroblastomaProcessed$feature.mat[,c("log2.n","log.hall")]
    n.folds <- 3
    uniq.folds <- 1:n.folds
    fold.vec <- sample(rep(uniq.folds, l=nrow(signal.features)))
    test.fold <- 1
    index.list <- list(
      train=fold.vec!=test.fold,
      test=fold.vec==test.fold)
    n.noise <- 40
    set.seed(1)
    noise.features <- matrix(
      rnorm(n.noise*nrow(signal.features)),
      nrow(signal.features), n.noise)
    input.features <- cbind(-54, 1, signal.features, noise.features)
    set.data <- lapply(index.list, function(is.set){
      L <- list(features=input.features[is.set,])
      for(denominator in c("count","rate")){
        L[[paste0(denominator, ".diffs")]] <- aum::aum_diffs_penalty(
          nb.err, rownames(L$features), denominator=denominator)
      }
      L
    })
    model <- with(set.data$train, aum::aum_linear_model_cv(
      features, count.diffs))
    plot(model)

    ## verify that the predictions are the same using either scaled or
    ## original features.
    with(set.data$train, rbind(
      predict.on.inputs=t(head(
        predict(model, features))),
      scale.then.predict=t(head(
        scale(features)[,model$keep] %*% model$weight.vec +
          model$intercept))))

    ## verify that the learned intercept results in min errors, at thresh=0.
    train.list <- with(set.data$train, aum::aum(
      count.diffs, predict(model, features)))
    plot(fp_before+fn_before ~ thresh, train.list$total_error)

    ## use rates instead of counts for computing AUM.
    set.seed(1)
    rate.model <- with(set.data$train, aum::aum_linear_model_cv(
      features, rate.diffs))
    plot(rate.model)

    ## alternative visualization including error bands and min loss.
    if(requireNamespace("ggplot2")){
      ggplot2::ggplot()+
        ggplot2::geom_ribbon(ggplot2::aes(
          step.number, ymin=aum_mean-aum_sd, ymax=aum_mean+aum_sd, fill=set),
          alpha=0.5,
          data=rate.model$set.loss)+
        ggplot2::geom_line(ggplot2::aes(
          step.number, aum_mean, color=set),
          data=rate.model$set.loss)+
        ggplot2::geom_point(ggplot2::aes(
          step.number, aum_mean, color=set),
          data=rate.model$set.loss[, .SD[which.min(aum_mean)], by=set])+
        ggplot2::scale_y_log10()
    }
    
    ## alternative visualization showing each fold.
    if(requireNamespace("ggplot2")){
      ggplot2::ggplot()+
        ggplot2::geom_line(ggplot2::aes(
          step.number, aum, color=set),
          data=rate.model$fold.loss)+
        ggplot2::geom_point(ggplot2::aes(
          step.number, aum, color=set),
          data=rate.model$fold.loss[
          , .SD[which.min(aum)], by=.(valid.fold, set)])+
        ggplot2::scale_y_log10()+
        ggplot2::facet_grid(. ~ valid.fold, labeller = "label_both")
    }

    ## compute test ROC curves, and compare against L1 regularized
    ## linear model with squared hinge loss.
    target.dt <- penaltyLearning::targetIntervals(nb.err, "example")
    target.mat <- target.dt[
      rownames(set.data$train$features),
      cbind(min.log.lambda, max.log.lambda),
      on="example"]
    named.features <- as.matrix(data.frame(input.features))
    ircv <- penaltyLearning::IntervalRegressionCV(
      named.features[rownames(set.data$train$features),], target.mat)
    pred.list <- with(set.data$test, list(
      IRCV=-predict(ircv, named.features[rownames(features),]),
      AUM.count=predict(model, features),
      AUM.rate=predict(rate.model, features),
      zero=rep(0, nrow(features))))
    roc.dt <- data.table(pred.name=names(pred.list))[, {
      aum::aum(set.data$test$rate.diffs, pred.list[[pred.name]])$total_error
    }, by=pred.name][, tp_before := 1-fn_before][]
    setkey(roc.dt, pred.name, thresh)
    if(requireNamespace("ggplot2")){
      pred.dt <- roc.dt[0 < thresh, .SD[1], by=pred.name]
      ggplot2::ggplot()+
        ggplot2::ggtitle("Test set ROC curves, dot for predicted threshold")+
        ggplot2::geom_path(ggplot2::aes(
          fp_before, tp_before, color=pred.name),
          data=roc.dt)+
        ggplot2::geom_point(ggplot2::aes(
          fp_before, tp_before, color=pred.name),
          shape=21,
          fill="white",
          data=pred.dt)+
        ggplot2::coord_equal(xlim=c(0,0.25), ylim=c(0.75,1))+
        ggplot2::xlab("False Positive Rate")+
        ggplot2::ylab("True Positive Rate")
    }
    ## first and last row of each pred.
    roc.dt[, .SD[c(1,.N)], by=pred.name]

    ## visualize area under min.
    roc.dt[, min_before := pmin(fp_before, fn_before)]
    roc.tall <- nc::capture_melt_single(
      roc.dt, 
      error.type="fp|fn|min",
      "_before",
      value.name="error.value")
    err.sizes <- c(
      fp=3,
      fn=2,
      min=1)
    err.colors <- c(
      fp="red",
      fn="deepskyblue",
      min="black")
    if(requireNamespace("ggplot2")){
      ggplot2::ggplot()+
        ggplot2::ggtitle("Train set AUM in green")+
        ggplot2::theme_bw()+
        ggplot2::geom_step(ggplot2::aes(
          thresh, error.value, color=error.type, linewidth=error.type),
          data=roc.tall)+
        ggplot2::geom_polygon(ggplot2::aes(
          thresh, min_before),
          fill="green",
          data=roc.dt)+
        ggplot2::facet_grid(pred.name ~ ., labeller="label_both")+
        ggplot2::xlab("Constant added to predicted values")+
        ggplot2::ylab("Error rate")+
        ggplot2::scale_color_manual(values=err.colors)+
        ggplot2::scale_size_manual(values=err.sizes)
    }

  }
  
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
  maxIterations=nrow(feature.list$subtrain)
### max number of iterations of exact line search, default is number
### of subtrain examples.
){
  weight.vec <- rep(0, ncol(feature.list$subtrain))
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
