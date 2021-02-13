#include <Rcpp.h>
#include <R.h>
#include "aum.h"
 
// [[Rcpp::export]]
Rcpp::List aum_interface
(const Rcpp::DataFrame error_df,
 const Rcpp::NumericVector pred_vec) {
  int pred_N = pred_vec.size();
  if(pred_N < 1){
    Rcpp::stop("need at least one prediction"); 
  }
  Rcpp::NumericVector out_aum(1);
  Rcpp::NumericMatrix out_deriv_mat(pred_N, 2);
  Rcpp::NumericVector max_pred = error_df["max_pred"];
  Rcpp::NumericVector fp = error_df["fp"];
  Rcpp::NumericVector fn = error_df["fn"];
  int err_N = error_df.nrow();
  Rcpp::NumericVector out_thresh(err_N);
  Rcpp::NumericVector out_fp_prb_diff(err_N);
  Rcpp::NumericVector out_fn_prb_diff(err_N);
  Rcpp::IntegerVector example = error_df["example"];
  int status = aum
    (&max_pred[0],
     &fp[0],
     &fn[0],
     &example[0],
     err_N,
     &pred_vec[0],
     pred_vec.size(),
     //inputs above, outputs below.
     &out_aum[0], &out_deriv_mat[0],
     &out_thresh[0],
     &out_fp_prb_diff[0],
     &out_fn_prb_diff[0]);     
  if(status != 0){
    Rcpp::stop("non-zero status"); 
  }
  return Rcpp::List::create
    (Rcpp::Named("aum", out_aum),
     Rcpp::Named("derivative_mat", out_deriv_mat)
     ) ;
}

