#define ERROR_FP_SHOULD_BE_ZERO_AT_NEGATIVE_INFINITY 1
#define ERROR_FN_SHOULD_BE_ZERO_AT_INFINITY 2

class Diffs{
  public:
  double fp_diff, fn_diff, fp_after, fn_before, min_before;
  Diffs();
};

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
   double *out_fn_prb_diff);
