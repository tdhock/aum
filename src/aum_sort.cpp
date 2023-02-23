#include "aum_sort.h"
#include <math.h>//isfinite
#include <algorithm>//std::sort

double get_min_thresh(const double *diff_vec, int err_N){
  double thresh = INFINITY;
  for(int err_i=0; err_i<err_N; err_i++){
    double abs_diff = abs(diff_vec[err_i]);
    if(abs_diff != 0 && abs_diff < thresh){
      thresh = abs_diff;
    }
  }
  return thresh/2;
}

// Main function for computing Area Under Minimum of False Positives
// and False Negatives. All pointer arguments must be arrays (of size
// indicated in comments below) that are allocated before calling
// aum_sort. Since we do std::sort on err_N elements the time
// complexity is O( err_N log err_N ).
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
  for(int err_i=0; err_i<err_N; err_i++){
    int example = err_example[err_i];
    if(example >= pred_N){
      return ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS;
    }
    if(example < 0){
      return ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_NON_NEGATIVE;
    }
    out_thresh[err_i] = err_pred[err_i] - pred_vec[example];
    out_indices[err_i] = err_i;
  }
  // Sort indices by threshold.
  std::sort
    (out_indices, out_indices+err_N,
     [&out_thresh](int left, int right){
       return out_thresh[left] < out_thresh[right];
     });
  double cumsum, cumsum_prev;
  const double *fp_or_fn_diff;
  double *out_this, *out_prev;
  double fp_or_fn_thresh;
  int first, err;
  double
    fp_diff_thresh=get_min_thresh(err_fp_diff, err_N),
    fn_diff_thresh=get_min_thresh(err_fn_diff, err_N);
  for(int err_type=0; err_type<2; err_type++){
    // Compute cumsums before AND after each threshold - fp starts at
    // zero and iterates forward, fn ends at zero and iterates backward.
    if(err_type == 0){//fn
      first = err_N - 1;
      sign = -1;
      fp_or_fn_thresh = fn_diff_thresh;
      fp_or_fn_diff = err_fn_diff;
      out_this = out_fn_before;
      out_prev = out_fn_after;
      err = ERROR_AUM_SORT_FN_SHOULD_BE_NON_NEGATIVE;
    }else{//fp
      first = 0;
      sign = 1;
      fp_or_fn_thresh = fp_diff_thresh;
      fp_or_fn_diff = err_fp_diff;
      out_this = out_fp_after;
      out_prev = out_fp_before;
      err = ERROR_AUM_SORT_FP_SHOULD_BE_NON_NEGATIVE;
    }
    cumsum = 0.0;
    cumsum_prev = 0.0;
    int i_prev = 0;
    for(int i=0; i<err_N; i++){
      int rank_i = first + sign*i;
      int err_i = out_indices[rank_i];
      cumsum += sign * fp_or_fn_diff[err_i];
      if(cumsum < -fp_or_fn_thresh){
	return err;
      }
      bool write_out = true;
      if(i == err_N-1){
        write_out = true;
      }else{
        int err_before = out_indices[rank_i+sign];
        write_out = out_thresh[err_i] != out_thresh[err_before];
      }
      if(write_out){
	for(int j=i_prev; j <= i; j++){
	  int rank_j = first + sign*j;
	  int err_j = out_indices[rank_j];
	  out_this[err_j] = cumsum;
	  out_prev[err_j] = cumsum_prev;
	}
	cumsum_prev = cumsum;
	i_prev = i+1;
      }
    }
  }
  // Compute AUM.
  for(int rank_i=1; rank_i<err_N; rank_i++){
    int err_i = out_indices[rank_i];
    int err_before = out_indices[rank_i-1];
    fp = out_fp_before[err_i];
    fn = out_fn_before[err_i];
    if(fp < fn){
      min = fp;
    }else{
      min = fn;
    }
    *out_aum += min*(out_thresh[err_i] - out_thresh[err_before]);
  }
  // Compute directional derivatives.
  for(int err_i=0; err_i<err_N; err_i++){
    int example = err_example[err_i];
    for(int out_col=0; out_col<2; out_col++){
      if(out_col==0){
	//for out_col=0, sign=-1, we join using out_thresh = min_thresh,
	//which means we need to take fp/fn/min from after the iterator.
	sign = -1;
	fp = out_fp_after[err_i];
	fn = out_fn_after[err_i];
      }else{//out_col=1, take before iterator.
	sign = 1;
	fp = out_fp_before[err_i];
	fn = out_fn_before[err_i];
      }
      if(fp < fn){
	min = fp;
      }else{
	min = fn;
      }
      fp_adj = fp + sign*err_fp_diff[err_i];
      fn_adj = fn + sign*err_fn_diff[err_i];
      if(fp_adj < fn_adj){
	min_adj = fp_adj;
      }else{
	min_adj = fn_adj;
      }
      out_deriv_mat[example+out_col*pred_N] += sign*(min_adj - min);
    }
  }
  return 0;
}
