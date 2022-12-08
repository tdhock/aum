\name{aum_linear_model_cv}
\alias{aum_linear_model_cv}
\title{aum linear model cv}
\description{Cross-validation for learning number of early stopping gradient
descent steps with exact line search, in linear model for
minimizing AUM.}
\usage{aum_linear_model_cv(feature.mat, 
    diff.dt, maxIterations = nrow(feature.mat), 
    improvement.thresh = NULL, 
    n.folds = 3)}
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
}

\value{Model trained with best number of iterations, represented as a
list of class aum_linear_model_cv with named elements: keep is a
logical vector telling which features should be kept before doing
matrix multiply of learned weight vector, weight.orig/weight.vec
and intercept.orig/intercept are the learned weights/intercepts
for the original/scaled feature space, fold.loss/set.loss are data
tables of loss values for the subtrain/validation sets, used for
selecting the best number of gradient descent steps.}

\author{Toby Hocking <toby.hocking@r-project.org> [aut, cre], Jadon Fowler [aut] (Contributed exact line search C++ code)}




\examples{

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

## verify that the predictions are the same using either scaled or
## original features.
rbind(
  predict.on.inputs=t(head(
    predict(model, input.features))),
  scale.then.predict=t(head(
    scale(input.features)[,model$keep] \%*\% model$weight.vec +
      model$intercept)))

## verify that the learned intercept results in min errors, at thresh=0.
train.list <- aum::aum(nb.diffs, predict(model, input.features))
plot(fp_before+fn_before ~ thresh, train.list$total_error)

## use rates instead of counts for computing AUM.
rate.diffs <- aum::aum_diffs_penalty(nb.err, rownames(input.features), denominator = "rate")
set.seed(1)
rate.model <- aum::aum_linear_model_cv(input.features, rate.diffs)
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

## compute ROC curves.
pred.list <- list(
  aum.count=predict(model, input.features),
  aum.rate=predict(rate.model, input.features),
  zero=rep(0, nrow(input.features)))
roc.dt <- data.table(pred.name=names(pred.list))[, {
  aum::aum(rate.diffs, pred.list[[pred.name]])$total_error
}, by=pred.name][, tp_before := 1-fn_before][]
setkey(roc.dt, pred.name, thresh)
if(requireNamespace("ggplot2")){
  pred.dt <- roc.dt[0 < thresh, .SD[1], by=pred.name]
  ggplot2::ggplot()+
    ggplot2::ggtitle("Train set ROC curves, dot for predicted threshold")+
    ggplot2::geom_path(ggplot2::aes(
      fp_before, tp_before, color=pred.name),
      data=roc.dt)+
    ggplot2::geom_point(ggplot2::aes(
      fp_before, tp_before, color=pred.name),
      shape=21,
      fill="white",
      data=pred.dt)+
    ggplot2::coord_equal()+
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
      thresh, error.value, color=error.type, size=error.type),
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
