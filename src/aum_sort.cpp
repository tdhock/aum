#include "aum_sort.h"
#include <math.h>//isfinite
#include <stdlib.h>//qsort_r

int compare_indices(const void *left, const void *right, void *out_ptr){
  double *out_thresh = (double*) out_ptr;
  return out_thresh[* (int*)left] > out_thresh[* (int*)right];
}

// Main function for computing Area Under Minimum of False Positives
// and False Negatives. All pointer arguments must be arrays (of size
// indicated in comments below) that are allocated before calling
// aum. Since we do qsort_r on err_N elements the time complexity is
// O( err_N log err_N ). A C compiler can be used with this code file
// but the "cpp" suffix is for easy compilation with the other C++
// code in the R package. This code is slightly faster (by constant
// factors) than aum_map, but this code is also slightly more
// complicated (caller must pre-allocate all temporary arrays, etc).
int aum_sort
(const double *err_pred, //err_N
 const double *err_fp_diff, 
 const double *err_fn_diff, 
 const int *err_example, 
 const int err_N,
 const double *pred_vec, //pred_N
 const int pred_N,
 //inputs above, outputs below.
 int *out_indices, //err_N
 double *out_thresh,
 double *out_fp_before,
 double *out_fp_after,
 double *out_fn_before,
 double *out_fn_after,
 double *out_aum, //1
 double *out_deriv_mat //pred_N*2
 ){
  double fp, fn, min, fp_adj, fn_adj, min_adj;
  int sign;
  *out_aum = 0.0;
  for(int out_i=0; out_i<pred_N*2; out_i++){
    out_deriv_mat[out_i] = 0.0;
  }
  for(int pred_i=0; pred_i<pred_N; pred_i++){
    if(!isfinite(pred_vec[pred_i])){
      return ERROR_AUM_SORT_ALL_PREDICTIONS_SHOULD_BE_FINITE;
    }
  }
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    if(row_example >= pred_N){
      return ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS;
    }
    if(row_example < 0){
      return ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_NON_NEGATIVE;
    }
    double pred_val = pred_vec[row_example];
    out_thresh[row] = err_pred[row] - pred_val;
    out_indices[row] = row;
  }
  // Sort indices by threshold.
  qsort_r(out_indices, err_N, sizeof(int), compare_indices, out_thresh);
  double cumsum, cumsum_prev;
  const double *fp_or_fn_diff;
  double *out_this, *out_prev;
  int first, err;
  for(int err_type=0; err_type<2; err_type++){
    // Compute cumsums before AND after each threshold - fp starts at
    // zero and iterates forward, fn ends at zero and iterates backward.
    if(err_type == 0){//fn
      first = err_N - 1;
      sign = -1;
      fp_or_fn_diff = err_fn_diff;
      out_this = &out_fn_before[0];
      out_prev = &out_fn_after[0];
      err = ERROR_AUM_SORT_FN_SHOULD_BE_NON_NEGATIVE;
    }else{//fp
      first = 0;
      sign = 1;
      fp_or_fn_diff = err_fp_diff;
      out_this = &out_fp_after[0];
      out_prev = &out_fp_before[0];
      err = ERROR_AUM_SORT_FP_SHOULD_BE_NON_NEGATIVE;
    }
    cumsum = 0.0;
    cumsum_prev = 0.0;
    int k_prev = 0;
    for(int k=0; k<err_N; k++){
      int i = first + sign*k;
      int row_i = out_indices[i];
      cumsum += sign * fp_or_fn_diff[row_i];
      if(cumsum < 0){
	return err;
      }
      if(k == err_N-1 || out_thresh[row_i] != out_thresh[out_indices[i+sign]]){
	for(int j=k_prev; j <= k; j++){
	  int row_j = out_indices[first+sign*j];
	  out_this[row_j] = cumsum;
	  out_prev[row_j] = cumsum_prev;
	}
	cumsum_prev = cumsum;
	k_prev = k+1;
      }
    }
  }
  // Compute AUM.
  for(int i=1; i<err_N; i++){
    int row_i = out_indices[i];
    int row_before = out_indices[i-1];
    fp = out_fp_before[row_i];
    fn = out_fn_before[row_i];
    if(fp < fn){
      min = fp;
    }else{
      min = fn;
    }
    *out_aum += min*(out_thresh[row_i] - out_thresh[row_before]);
  }
  // Compute directional derivatives.
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    for(int out_col=0; out_col<2; out_col++){
      if(out_col==0){
	//for out_col=0, sign=-1, we join using out_thresh = min_thresh,
	//which means we need to take fp/fn/min from after the iterator.
	sign = -1;
	fp = out_fp_after[row];
	fn = out_fn_after[row];
      }else{//out_col=1, take before iterator.
	sign = 1;
	fp = out_fp_before[row];
	fn = out_fn_before[row];
      }
      if(fp < fn){
	min = fp;
      }else{
	min = fn;
      }
      fp_adj = fp + sign*err_fp_diff[row];
      fn_adj = fn + sign*err_fn_diff[row];
      if(fp_adj < fn_adj){
	min_adj = fp_adj;
      }else{
	min_adj = fn_adj;
      }
      out_deriv_mat[row_example+out_col*pred_N] += sign*(min_adj - min);
    }
  }
  return 0;
}
