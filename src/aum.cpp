#include "aum.h"
#include <map>

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
int aum
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
  TotalsMap totals_map;
  Totals zero_diffs;
  for(int row=0; row<err_N; row++){
    int row_example = err_example[row];
    out_thresh[row] = err_pred[row] - pred_vec[row_example];
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
