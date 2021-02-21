#define ERROR_AUM_SORT_FP_SHOULD_BE_NON_NEGATIVE 1
#define ERROR_AUM_SORT_FN_SHOULD_BE_NON_NEGATIVE 2
#define ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS 3
#define ERROR_AUM_SORT_EXAMPLE_SHOULD_BE_NON_NEGATIVE 4
#define ERROR_AUM_SORT_ALL_PREDICTIONS_SHOULD_BE_FINITE 5

int aum_sort
  (const double *err_pred,
   const double *err_fp,
   const double *err_fn,
   const int *err_example,
   const int err_N,
   const double *pred_vec,
   const int pred_N,
   //inputs above, outputs below.
   int *out_indices, //err_N
   double *out_thresh,
   double *out_fp_before,
   double *out_fp_after,
   double *out_fn_before,
   double *out_fn_after,
   double *out_aum,
   double *out_deriv_mat
   );
