#include <Rcpp.h>
#include "aum_sort.h"
#include "aum_map.h"
#include <omp.h> 

// [[Rcpp::export]]
Rcpp::NumericVector stl_sort(Rcpp::NumericVector x) {
  Rcpp::NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}
 
// [[Rcpp::export]]
Rcpp::NumericVector stl_set_insert(Rcpp::NumericVector x) {
  std::set<double> s;
  for(int i=0; i<x.size(); i++){
    s.insert(x[i]);
  }
  std::set<double>::iterator it = s.begin();
  Rcpp::NumericVector y = clone(x);
  for(int i=0; i<x.size(); i++){
    y[i] = *it;
    it++;
  }
  return y;
}
 
// [[Rcpp::export]]
Rcpp::List aum_map_interface
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
  Rcpp::NumericVector out_aum(1);
  Rcpp::NumericMatrix out_deriv_mat(pred_N, 2);
  int status = aum_map
    (&err_pred[0],
     &err_fp_diff[0],
     &err_fn_diff[0],
     &err_example[0],
     err_N,
     &pred_vec[0],
     pred_vec.size(),
     //inputs above, outputs below.
     &out_aum[0], &out_deriv_mat[0]);
  if(status == ERROR_AUM_MAP_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS){
    Rcpp::stop("example should be less than number of predictions"); 
  }
  if(status == ERROR_AUM_MAP_EXAMPLE_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("example should be non-negative");
  }
  if(status == ERROR_AUM_MAP_FN_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fn should be non-negative"); 
  }
  if(status == ERROR_AUM_MAP_FP_SHOULD_BE_NON_NEGATIVE){
    Rcpp::stop("fp should be non-negative"); 
  }
  if(status == ERROR_AUM_MAP_ALL_PREDICTIONS_SHOULD_BE_FINITE){
    Rcpp::stop("all predictions should be finite");
  }
  return Rcpp::List::create
    (Rcpp::Named("aum", out_aum),
     Rcpp::Named("derivative_mat", out_deriv_mat)
     ) ;
}

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
     Rcpp::Named("derivative_mat", out_deriv_mat)
     ) ;
}

// [[Rcpp::export]]
void multithread
(const Rcpp::DataFrame err_df,
 const Rcpp::NumericVector pred_vec,
 int threads){
#pragma omp parallel for
  for(int i=0; i<threads; i++){
    aum_sort_interface(err_df, pred_vec);
  }
}
