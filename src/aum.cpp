#include "aum.h"
#include <math.h>//isfinite
#include <map>
#include <algorithm>//std::sort

Totals::Totals(){
  fp_diff = 0;
  fn_diff = 0;
}

// Need to use a map rather than a set because map values are mutable
// whereas set elements are not.
typedef std::map<double, Totals> TotalsMap;

// Main function for computing Area Under Minimum of False Positives
// and False Negatives. All pointer arguments must be arrays (of size
// indicated in comments below) that are allocated before calling
// aum. Since err_N >= pred_N the complexity is O( err_N log err_N ),
// because we use log-time STL map methods (insert, find) err_N times.
int aum_map
(const double *err_pred, //err_N
 const double *err_fp_diff, //err_N
 const double *err_fn_diff, //err_N
 const int *err_example, //err_N
 const int err_N,
 const double *pred_vec, //pred_N
 const int pred_N,
 //inputs above, outputs below.
 double *out_thresh, //err_N
 double *out_aum, //1
 double *out_deriv_mat //pred_N*2
 ){
  *out_aum = 0.0;
  for(int out_i=0; out_i<pred_N*2; out_i++){
    out_deriv_mat[out_i] = 0.0;
  }
  for(int pred_i=0; pred_i<pred_N; pred_i++){
    if(!isfinite(pred_vec[pred_i])){
      return ERROR_ALL_PREDICTIONS_SHOULD_BE_FINITE;
    }
  }
  TotalsMap totals_map;
  Totals zero_diffs;
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    if(row_example >= pred_N){
      return ERROR_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS;
    }
    if(row_example < 0){
      return ERROR_EXAMPLE_SHOULD_BE_NON_NEGATIVE;
    }
    double pred_val = pred_vec[row_example];
    out_thresh[row] = err_pred[row] - pred_val;
    std::pair<double,Totals> to_insert(out_thresh[row],zero_diffs);
    std::pair<TotalsMap::iterator,bool> ret;
    ret = totals_map.insert(to_insert);
    ret.first->second.fp_diff += err_fp_diff[row];
    ret.first->second.fn_diff += err_fn_diff[row];
  }
  double cumsum;
  cumsum = 0.0;
  for(TotalsMap::iterator it=totals_map.begin();
      it != totals_map.end(); it++){
    cumsum += it->second.fp_diff;
    if(cumsum < 0){
      return ERROR_FP_SHOULD_BE_NON_NEGATIVE;
    }
    it->second.fp_after = cumsum;
  }
  cumsum = 0.0;
  for(TotalsMap::reverse_iterator it=totals_map.rbegin();
      it != totals_map.rend(); it++){
    cumsum -= it->second.fn_diff;
    if(cumsum < 0){
      return ERROR_FN_SHOULD_BE_NON_NEGATIVE;
    }
    it->second.fn_before = cumsum;
    TotalsMap::reverse_iterator before = next(it);
    double min;
    if(before == totals_map.rend()){
      min = 0;
    }else{
      double fp_before = before->second.fp_after;
      if(cumsum < fp_before){
	min = cumsum;
      }else{
	min = fp_before;
      }
      *out_aum += min*(it->first - before->first);
    }
    it->second.min_before = min;
  }
  double sign;
  for(int out_col=0; out_col<2; out_col++){
    if(out_col==0){
      sign = -1;
    }else{
      sign = 1;
    }
    for(int row=0; row<err_N; row++){
      int row_example = err_example[row];
      TotalsMap::iterator it = totals_map.find(out_thresh[row]);
      double fp, fn, min, fp_adj, fn_adj, min_adj;
      if(out_col==0){
	//for out_col=0, sign=-1, we join using out_thresh = min_thresh,
	//which means we need to take fp/fn/min from after the iterator.
	fp = it->second.fp_after;
	if(next(it) == totals_map.end()){
	  fn = 0;
	  min = 0;
	}else{
	  fn = next(it)->second.fn_before;
	  min = next(it)->second.min_before;
	}
      }else{//out_col=1, take before iterator.
	if(it == totals_map.begin()){
	  fp = 0;
	}else{
	  fp = prev(it)->second.fp_after;
	}
	fn = it->second.fn_before;
	min = it->second.min_before;
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

int aum_sort
(const double *err_pred, //err_N
 const double *err_fp_diff, //err_N
 const double *err_fn_diff, //err_N
 const int *err_example, //err_N
 const int err_N,
 const double *pred_vec, //pred_N
 const int pred_N,
 //inputs above, outputs below.
 int *out_indices,
 double *out_thresh, //err_N
 double *out_fp_before,
 double *out_fp_after,
 double *out_fn_before,
 double *out_fn_after,
 double *out_aum, //1
 double *out_deriv_mat //pred_N*2
 ){
  *out_aum = 0.0;
  for(int out_i=0; out_i<pred_N*2; out_i++){
    out_deriv_mat[out_i] = 0.0;
  }
  for(int pred_i=0; pred_i<pred_N; pred_i++){
    if(!isfinite(pred_vec[pred_i])){
      return ERROR_ALL_PREDICTIONS_SHOULD_BE_FINITE;
    }
  }
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    if(row_example >= pred_N){
      return ERROR_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS;
    }
    if(row_example < 0){
      return ERROR_EXAMPLE_SHOULD_BE_NON_NEGATIVE;
    }
    double pred_val = pred_vec[row_example];
    out_thresh[row] = err_pred[row] - pred_val;
    out_indices[row] = row;
  }
  std::sort
    (out_indices,
     out_indices+err_N,
     [&out_thresh](int left, int right){
       return out_thresh[left] < out_thresh[right];
    }
  );
  double cumsum_after, cumsum_before;
  cumsum_after = 0.0;
  cumsum_before = 0.0;
  int i_prev = 0;
  for(int i=0; i<err_N; i++){
    int row_i = out_indices[i];
    cumsum_after += err_fp_diff[row_i];
    if(cumsum_after < 0){
      return ERROR_FP_SHOULD_BE_NON_NEGATIVE;
    }
    int i_next = i+1;
    int row_next = out_indices[i_next];
    if(i == err_N-1 || out_thresh[row_i] != out_thresh[row_next]){
      for(int j=i_prev; j <= i; j++){
	int row_j = out_indices[j];
	out_fp_after[row_j] = cumsum_after;
	out_fp_before[row_j] = cumsum_before;
      }
      cumsum_before = cumsum_after;
      i_prev = i+1;
    }
  }
  cumsum_after = 0.0;
  cumsum_before = 0.0;
  i_prev = err_N-1;
  for(int i=err_N-1; i >= 0; i--){
    int row_i = out_indices[i];
    cumsum_before -= err_fn_diff[row_i];
    if(cumsum_before < 0){
      return ERROR_FN_SHOULD_BE_NON_NEGATIVE;
    }
    int i_next=i-1;
    int row_next = out_indices[i_next];
    if(i == 0 || out_thresh[row_next] != out_thresh[row_i]){
      for(int j=i; j <= i_prev; j++){
	int row_j = out_indices[j];
	out_fn_before[row_j] = cumsum_before;
	out_fn_after[row_j] = cumsum_after;
      }
      cumsum_after = cumsum_before;
      i_prev = i-1;
    }
  }
  double fp, fn, min, fp_adj, fn_adj, min_adj, sign;
  for(int i=1; i<err_N; i++){
    int row_i = out_indices[i];
    int row_before = out_indices[i-1];
    if(out_fp_before[row_i] < out_fn_before[row_i]){
      min = out_fp_before[row_i];
    }else{
      min = out_fn_before[row_i];
    }
    *out_aum += min*(out_thresh[row_i] - out_thresh[row_before]);
  }
  for(int out_col=0; out_col<2; out_col++){
    if(out_col==0){
      sign = -1;
    }else{
      sign = 1;
    }
    for(int row=0; row<err_N; row++){
      int row_example = err_example[row];
      if(out_col==0){
	//for out_col=0, sign=-1, we join using out_thresh = min_thresh,
	//which means we need to take fp/fn/min from after the iterator.
	fp = out_fp_after[row];
	fn = out_fn_after[row];
      }else{//out_col=1, take before iterator.
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
