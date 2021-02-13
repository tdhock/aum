#include "aum.h"
#include "math.h"
#include <map>
#include <stdio.h>

Diffs::Diffs(){
  fp_diff = 0;
  fn_diff = 0;
}

typedef std::map<double, Diffs> DiffsMap;
 
int aum
  (const double *err_max_pred,
   const double *err_fp,
   const double *err_fn,
   const int *err_example,
   const int err_N,
   const double *pred_vec,
   const int pred_N,
   //inputs above, outputs below.
   double *out_aum,
   double *out_deriv_mat,
   double *out_thresh,
   double *out_fp_prb_diff,
   double *out_fn_prb_diff){
  DiffsMap thresh_map;
  double cumsum;
  *out_aum = 0.0;
  Diffs zero_diffs;
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    out_thresh[row] = err_max_pred[row] - pred_vec[row_example];
    bool last_row = row==err_N-1;
    bool another_example_next =
      (!last_row) && (err_example[row] < err_example[row+1]);
    bool last_row_in_example = last_row || another_example_next;
    bool first_row = row==0;
    bool another_example_before =
      (!first_row) && (err_example[row-1] < err_example[row]);
    bool first_row_in_example = first_row || another_example_before;
    if(first_row_in_example && err_fp[row] != 0){
      return ERROR_FP_SHOULD_BE_ZERO_AT_NEGATIVE_INFINITY;
    }
    if(last_row_in_example){
      if(err_fn[row] != 0){
	return ERROR_FN_SHOULD_BE_ZERO_AT_INFINITY;
      }
      out_fp_prb_diff[row] = NAN;
      out_fn_prb_diff[row] = NAN;
    }else{
      out_fp_prb_diff[row] = err_fp[row+1]-err_fp[row];
      out_fn_prb_diff[row] = err_fn[row+1]-err_fn[row];
      std::pair<double,Diffs> to_insert(out_thresh[row],zero_diffs);
      std::pair<DiffsMap::iterator,bool> ret;
      ret = thresh_map.insert(to_insert);
      ret.first->second.fp_diff += out_fp_prb_diff[row];
      ret.first->second.fn_diff += out_fn_prb_diff[row];
    }
  }
  cumsum = 0.0;
  for(DiffsMap::iterator it=thresh_map.begin(); it != thresh_map.end(); it++){
    cumsum += it->second.fp_diff;
    it->second.fp_after = cumsum;
  }
  cumsum = 0.0;
  for(DiffsMap::reverse_iterator it=thresh_map.rbegin(); it != thresh_map.rend(); it++){
    cumsum -= it->second.fn_diff;
    it->second.fn_before = cumsum;
    DiffsMap::reverse_iterator before = next(it);
    if(before != thresh_map.rend()){
      double fp_before = before->second.fp_after;
      double min;
      if(cumsum < fp_before){
	min = cumsum;
      }else{
	min = fp_before;
      }
      it->second.min_before = min;
      *out_aum += min*(it->first - before->first);
    }
  }
  for(DiffsMap::iterator it=thresh_map.begin(); it != thresh_map.end(); it++){
    printf("%10.2f fp_diff=%3.0f fn_diff=%3.0f fn_before=%3.0f fp_after=%3.0f\n", it->first, it->second.fp_diff, it->second.fn_diff, it->second.fn_before, it->second.fp_after);
  }
  return 0;
}
