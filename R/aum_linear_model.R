aum_linear_model_cv <- structure(function(feature.mat, diff.dt, maxIterations=nrow(feature.mat), improvement.thresh=0.1, n.folds=3){
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
### Model trained with best number of iterations.
}, ex=function(){

  ## learn a model for a real changepoint data set.
  data(neuroblastomaProcessed, package="penaltyLearning", envir=environment())
  nb.err <- with(neuroblastomaProcessed$errors, data.frame(
    example=paste0(profile.id, ".", chromosome),
    min.lambda,
    max.lambda,
    fp, fn))
  signal.features <- neuroblastomaProcessed$feature.mat[,c("log2.n","log.hall")]
  n.noise <- 40
  set.seed(1)
  noise.features <- matrix(
    rnorm(n.noise*nrow(signal.features)),
    nrow(signal.features), n.noise)
  input.features <- cbind(-54, 1, signal.features, noise.features)
  nb.diffs <- aum::aum_diffs_penalty(nb.err, rownames(input.features))
  model <- aum::aum_linear_model_cv(input.features, nb.diffs)
  plot(model)

  ## alternative visualization including error bands and min loss.
  if(requireNamespace("ggplot2")){
    ggplot2::ggplot()+
      ggplot2::geom_ribbon(ggplot2::aes(
        step.number, ymin=aum_mean-aum_sd, ymax=aum_mean+aum_sd, fill=set),
        alpha=0.5,
        data=model$set.loss)+
      ggplot2::geom_line(ggplot2::aes(
        step.number, aum_mean, color=set),
        data=model$set.loss)+
      ggplot2::geom_point(ggplot2::aes(
        step.number, aum_mean, color=set),
        data=model$set.loss[, .SD[which.min(aum_mean)], by=set])+
      ggplot2::scale_y_log10()
  }

  ## verify that the predictions are the same using either scaled or
  ## original features.
  rbind(
    predict.on.inputs=t(head(
      predict(model, input.features))),
    scale.then.predict=t(head(
      scale(input.features)[,model$keep] %*% model$weight.vec +
        model$intercept)))

  ## verify that the learned intercept results in min errors, at thresh=0.
  train.list <- aum::aum(nb.diffs, predict(model, input.features))
  plot(fp_before+fn_before ~ thresh, train.list$total_error)
  
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

aum_linear_model <- function(feature.list, diff.list, max.steps=NULL, improvement.thresh=NULL, maxIterations=nrow(feature.list$subtrain)){
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
}  

plot.aum_linear_model <- function(x, ...){
  lattice::xyplot(
    aum ~ step.number, x$loss, 
    groups=set, type="l", 
    auto.key=list(space="right", points=FALSE, lines=TRUE))
}
