#include <Rcpp.h>
#include "aum.h"
 
// [[Rcpp::export]]
Rcpp::List aum_interface
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
  Rcpp::NumericVector out_thresh(err_N);
  Rcpp::NumericVector out_aum(1);
  Rcpp::NumericMatrix out_deriv_mat(pred_N, 2);
  int status = aum
    (&err_pred[0],
     &err_fp_diff[0],
     &err_fn_diff[0],
     &err_example[0],
     err_N,
     &pred_vec[0],
     pred_vec.size(),
     //inputs above, outputs below.
     &out_thresh[0],
     &out_aum[0], &out_deriv_mat[0]);
  if(status == ERROR_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS){
    Rcpp::stop("example should be less than number of predictions"); 
  }
  if(status == ERROR_EXAMPLE_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("example should be non-negative");
  }
  if(status == ERROR_FN_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fn should be non-negative"); 
  }
  if(status == ERROR_FP_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fp should be non-negative"); 
  }
  if(status != 0){
    Rcpp::stop("non-zero status"); 
  }
  return Rcpp::List::create
    (Rcpp::Named("aum", out_aum),
     Rcpp::Named("derivative_mat", out_deriv_mat)
     ) ;
}

