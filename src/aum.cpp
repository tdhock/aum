#include "aum.h"
#include "math.h"
#include <map>
#include <stdio.h>

Totals::Totals(){
  fp_diff = 0;
  fn_diff = 0;
}

// Need to use a map rather than a set because map values are mutable
// whereas set elements are not.
typedef std::map<double, Totals> TotalsMap;
 
int aum
  (const double *err_pred,
   const double *err_fp_diff,
   const double *err_fn_diff,
   const int *err_example,
   const int err_N,
   const double *pred_vec,
   const int pred_N,
   //inputs above, outputs below.
   double *out_thresh,
   double *out_aum,
   double *out_deriv_mat){
  TotalsMap totals_map;
  double cumsum;
  *out_aum = 0.0;
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
  cumsum = 0.0;
  for(TotalsMap::iterator it=totals_map.begin(); it != totals_map.end(); it++){
    cumsum += it->second.fp_diff;
    if(cumsum < 0){
      return ERROR_FP_SHOULD_BE_NON_NEGATIVE;
    }
    it->second.fp_after = cumsum;
  }
  cumsum = 0.0;
  for(TotalsMap::reverse_iterator it=totals_map.rbegin(); it != totals_map.rend(); it++){
    cumsum -= it->second.fn_diff;
    if(cumsum < 0){
      return ERROR_FN_SHOULD_BE_NON_NEGATIVE;
    }
    it->second.fn_before = cumsum;
    TotalsMap::reverse_iterator before = next(it);
    if(before != totals_map.rend()){
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
  for(TotalsMap::iterator it=totals_map.begin(); it != totals_map.end(); it++){
    printf("%10.2f fp_diff=%3.0f fn_diff=%3.0f fn_before=%3.0f fp_after=%3.0f\n", it->first, it->second.fp_diff, it->second.fn_diff, it->second.fn_before, it->second.fp_after);
  }
  return 0;
}
