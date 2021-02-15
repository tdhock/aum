aum <- structure(function
### Compute the Area Under Minimum of False Positives and False
### Negatives, and its directional derivatives.
(error.diff.df,
### data.frame of error differences. There should be one row for each
### change in error functions. "example" column indicates example ID
### from 1 to N, "pred" column indicates predicted value where there
### is a change in the error function(s), "fp_diff" and "fn_diff"
### columns indicate differences in false positives and false
### negatives at that predicted value. Note that this representation
### assumes that each error function has fp=0 at pred=-Inf and fn=0 at
### pred=Inf.
  pred.vec
### numeric vector of N predicted values.
){
  if(is.integer(error.diff.df$example)){
    error.diff.df$example <- error.diff.df$example-1L
  }else{
    stop("error.diff.df$example is non-integer which is not yet supported")
  }
  aum_interface(error.diff.df, pred.vec)
### Named list of two items: aum is numeric scalar loss value,
### derivative_mat is N x 2 matrix of directional derivatives (first
### column is derivative from left, second column is derivative from
### right).
}, ex=function(){

  (bin.diffs <- aum::aum_diffs_binary(c(0,1)))
  aum::aum(bin.diffs, c(-1,1))
  aum::aum(bin.diffs, c(0,0))
  aum::aum(bin.diffs, c(1,-1))
  
})
