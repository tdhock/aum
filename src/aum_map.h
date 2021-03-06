#define ERROR_AUM_MAP_FP_SHOULD_BE_NON_NEGATIVE 1
#define ERROR_AUM_MAP_FN_SHOULD_BE_NON_NEGATIVE 2
#define ERROR_AUM_MAP_EXAMPLE_SHOULD_BE_LESS_THAN_NUMBER_OF_PREDICTIONS 3
#define ERROR_AUM_MAP_EXAMPLE_SHOULD_BE_NON_NEGATIVE 4
#define ERROR_AUM_MAP_ALL_PREDICTIONS_SHOULD_BE_FINITE 5

class Totals{
  public:
  double fp_diff, fn_diff, fp_after, fn_before, min_before;
  Totals();
};

int aum_map
  (const double *err_pred,
   const double *err_fp,
   const double *err_fn,
   const int *err_example,
   const int err_N,
   const double *pred_vec,
   const int pred_N,
   //inputs above, outputs below.
   double *out_aum,
   double *out_deriv_mat
   );

