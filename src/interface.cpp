#include <Rcpp.h>
#include "aum_sort.h"
#include "aum_line_search.h"

// [[Rcpp::export]]
Rcpp::List aum_sort_interface
(const Rcpp::DataFrame err_df,
 const Rcpp::NumericVector pred_vec) {
  int pred_N = pred_vec.size();
  if(pred_N < 1){
    Rcpp::stop("need at least one prediction"); 
  }
  Rcpp::NumericVector err_pred = err_df["pred"];
  Rcpp::NumericVector err_fp_diff = err_df["fp_diff"];
  Rcpp::NumericVector err_fn_diff = err_df["fn_diff"];
  Rcpp::IntegerVector err_example = err_df["example"];
  int err_N = err_df.nrow();
  Rcpp::IntegerVector out_indices(err_N);
  Rcpp::NumericVector out_thresh(err_N);
  Rcpp::NumericVector out_fp_before(err_N);
  Rcpp::NumericVector out_fp_after(err_N);
  Rcpp::NumericVector out_fn_before(err_N);
  Rcpp::NumericVector out_fn_after(err_N);
  Rcpp::NumericVector out_aum(1);
  Rcpp::NumericMatrix out_deriv_mat(pred_N, 2);
  int status = aum_sort
    (&err_pred[0],
     &err_fp_diff[0],
     &err_fn_diff[0],
     &err_example[0],
     err_N,
     &pred_vec[0],
     pred_vec.size(),
     //inputs above, outputs below.
     &out_indices[0],
     &out_thresh[0],
     &out_fp_before[0],
     &out_fp_after[0],
     &out_fn_before[0],
     &out_fn_after[0],
     &out_aum[0],
     &out_deriv_mat[0]);
  if(status == ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS){
    Rcpp::stop("example should be less than number of predictions"); 
  }
  if(status == ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("example should be non-negative");
  }
  if(status == ERROR_AUM_SORT_FN_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fn should be non-negative"); 
  }
  if(status == ERROR_AUM_SORT_FP_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fp should be non-negative"); 
  }
  if(status == ERROR_AUM_SORT_ALL_PREDICTIONS_SHOULD_BE_FINITE){
    Rcpp::stop("all predictions should be finite");
  }
  return Rcpp::List::create
    (Rcpp::Named("aum", out_aum),
     Rcpp::Named("derivative_mat", out_deriv_mat),
     Rcpp::Named
     ("total_error",
      Rcpp::DataFrame::create
      (Rcpp::Named("thresh", out_thresh),
       Rcpp::Named("fp_before", out_fp_before),
       Rcpp::Named("fn_before", out_fn_before))));
}

// [[Rcpp::export]]
Rcpp::DataFrame aumLineSearch(const Rcpp::DataFrame df, int maxIterations) {
    // extract columns from dataframe
    Rcpp::NumericVector fpDiff = df["fp.diff"];
    Rcpp::NumericVector fnDiff = df["fn.diff"];
    Rcpp::NumericVector intercept = df["intercept"];
    Rcpp::NumericVector slope = df["slope"];
    int lineCount = df.nrow();
    Rcpp::NumericVector stepSizeVec(maxIterations, -1.0);
    Rcpp::NumericVector aumVec(maxIterations, -1.0);
    Rcpp::NumericVector aumSlopeAfterStepVec(maxIterations, -100.0);
    Rcpp::NumericVector aucAtStepVec(maxIterations, -1.0);
    Rcpp::NumericVector aucAfterStepVec(maxIterations, -1.0);
    Rcpp::IntegerVector intersectionCountVec(maxIterations, -1);
    Rcpp::IntegerVector intervalCountVec(maxIterations, -1);
    Rcpp::IntegerVector qSizeVec(maxIterations, -1);
    int status = lineSearch(
            &intercept[0],
            &slope[0],
            lineCount,
            &fpDiff[0],
            &fnDiff[0],
            maxIterations,
            &stepSizeVec[0],
            &aumVec[0],
            &aumSlopeAfterStepVec[0],
            &aucAtStepVec[0],
            &aucAfterStepVec[0],
            &intersectionCountVec[0],
            &intervalCountVec[0],
	    &qSizeVec[0]
    );
    if(status == ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING){
      Rcpp::stop("intercepts should be non-decreasing");
    }
    if(status == ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS){
      Rcpp::stop("slopes should be increasing for equal intercepts");
    }
    if(status == ERROR_LINE_SEARCH_MAX_FP_SHOULD_BE_POSITIVE){
      Rcpp::stop("max FP should be positive");
    }
    if(status == ERROR_LINE_SEARCH_MAX_FN_SHOULD_BE_POSITIVE){
      Rcpp::stop("max FN should be positive");
    }
    return Rcpp::DataFrame::create
      (Rcpp::Named("step.size", stepSizeVec),
       Rcpp::Named("aum", aumVec), 
       Rcpp::Named("aum.slope.after", aumSlopeAfterStepVec), 
       Rcpp::Named("auc", aucAtStepVec),
       Rcpp::Named("auc.after", aucAfterStepVec),
       Rcpp::Named("intersections", intersectionCountVec),
       Rcpp::Named("intervals", intervalCountVec),
       Rcpp::Named("q.size", qSizeVec));
}
