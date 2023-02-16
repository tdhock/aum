#include "aum_line_search.h"
using namespace std;

bool Point::isFinite() const {
    return isfinite(x) && isfinite(y);
}

Point intersect
(double a_intercept, double b_intercept, 
 double a_slope, double b_slope) {
    if (a_slope == b_slope) return Point{INFINITY, INFINITY};
    double x = (b_intercept - a_intercept) / (a_slope - b_slope);
    double y = a_intercept + a_slope * x;
    return Point{x, y};
}

class Intervals {
  public:
  int high_id, high_rank, low_rank, n_intervals;
  vector<int> *rank_from_id;
  Intervals(int high_id_, int n_intervals_, vector<int> *rank_from_id_){
    high_id = high_id_;
    n_intervals = n_intervals_;
    rank_from_id = rank_from_id_;
  }
  void set_ranks(){
    high_rank = (*rank_from_id)[high_id];
    low_rank = high_rank - n_intervals;
  }
};

class ThreshIntervals {
  public:
  map<double,Intervals> thresh_intervals_map;
  vector<int> *id_from_rank, *rank_from_id;
  ThreshIntervals
  (double thresh, 
   int high_rank,
   vector<int> *id_from_rank_,
   vector<int> *rank_from_id_
   ){
    id_from_rank = id_from_rank_;
    rank_from_id = rank_from_id_;
    int high_id = (*id_from_rank)[high_rank];
    thresh_intervals_map.insert
      (pair<double,Intervals>
       (thresh, Intervals(high_id,1,rank_from_id)));
  }
  void set_intervals_ranks(){
    for//intersection points. (tie-breaking type 1)
      (auto intervals_it = thresh_intervals_map.begin(); 
       intervals_it != thresh_intervals_map.end();
       intervals_it++){
      intervals_it->second.set_ranks();
    }
  }
  void add_interval(double thresh, int new_high_rank){
    auto it_after_or_same = thresh_intervals_map.lower_bound(thresh);
    int new_high_id = (*id_from_rank)[new_high_rank];
    if(it_after_or_same != thresh_intervals_map.end()){
      if(it_after_or_same->first == thresh){//same
        int old_high_rank = (*rank_from_id)[it_after_or_same->second.high_id];
        int old_low_rank = old_high_rank-it_after_or_same->second.n_intervals;
        if(old_high_rank+1 == new_high_rank){
          it_after_or_same->second.high_id = new_high_id;
        }else if(old_low_rank!=new_high_rank){
          return;
          //printf("WARNING: in add_interval, trying to alter existing/old Intervals (%d,%d) but new_high_rank=%d is not adjacent, this only happens when same as existing, no need to add again.\n", old_low_rank, old_high_rank, new_high_rank);
        }
        it_after_or_same->second.n_intervals++;
        return;
      }
    }
    thresh_intervals_map.insert
      (it_after_or_same, 
       pair<double,Intervals>
       (thresh, Intervals(new_high_id,1,rank_from_id)));
  }
};

class TotalAUC {
  public:
  vector<double> *FP, *FN, *M;
  vector<int> *id_from_rank, *rank_from_id;
  const double *slope_from_id;
  double value, aum_slope, maxFP, maxFN;
  int lineCount;
  TotalAUC
  (vector<double> *FP_, 
   double maxFP_,
   vector<double> *FN_, 
   double maxFN_,
   vector<double> *M_, 
   vector<int> *id_from_rank_,
   vector<int> *rank_from_id_,
   const double *slope_from_id_,
   int lineCount_){
    FP = FP_;
    maxFP = maxFP_;
    FN = FN_;
    maxFN = maxFN_;
    M = M_;
    id_from_rank = id_from_rank_;
    rank_from_id = rank_from_id_;
    slope_from_id = slope_from_id_;
    lineCount = lineCount_;
    zero();
  }
  void zero(){
    value = 0;
    aum_slope = 0;
  }
  void update_all(){
    update(1, lineCount, 1.0);
  }
  void zero_update_all(){
    zero();
    update_all();
  }
  double tpr(int rank){
    if(rank==lineCount)return 1;
    return 1-(*FN)[rank]/maxFN;
  }
  double fpr(int rank){
    if(rank==lineCount)return 1;
    return (*FP)[rank]/maxFP;
  }
  double m(int rank){
    if(rank==lineCount)return 0;
    return (*M)[rank];
  }
  double get_auc(int leftRank, int rightRank){
    double FPR_diff = fpr(rightRank)-fpr(leftRank);
    double TPR_sum = tpr(rightRank)+tpr(leftRank);
    return FPR_diff*TPR_sum/2;
  }
  void update
  (int first_rightRank, 
   int last_rightRank,
   double sign){
    for
      (int rightRank=first_rightRank; 
       rightRank <= last_rightRank; 
       rightRank++){
      int leftRank=rightRank-1;
      double roc_change = sign*get_auc(leftRank, rightRank);
      value += roc_change;
    }
    for
      (int rank=first_rightRank-1; 
       rank < last_rightRank; 
       rank++){
      int id = (*id_from_rank)[rank];
      double rank_slope = slope_from_id[id];
      double min_diff = m(rank)-m(rank+1);
      double slope_update = sign*rank_slope*min_diff;
      aum_slope += slope_update;
    }
  }
  double action(ThreshIntervals *TI_ptr, double sign){
    double more_auc_at_step = 0;
    for//intersection points. (tie-breaking type 1)
      (auto TI_it = TI_ptr->thresh_intervals_map.begin(); 
       TI_it != TI_ptr->thresh_intervals_map.end();
       TI_it++){
      int low_rank, high_rank;
      high_rank=TI_it->second.high_rank+1;
      low_rank=TI_it->second.low_rank+1;
      update(low_rank, high_rank, sign);
      more_auc_at_step += get_auc(low_rank-1, high_rank);
    }
    return more_auc_at_step;
  }
};

class Queue {
  public:
  map<double,ThreshIntervals> step_ThreshIntervals_map;
  vector<int> *id_from_rank, *rank_from_id;
  const double *intercept_from_id, *slope_from_id;
  Queue
  (vector<int> *id_from_rank_, 
   vector<int> *rank_from_id_,
   const double *intercept_from_id_,
   const double *slope_from_id_){
    id_from_rank = id_from_rank_;
    rank_from_id = rank_from_id_;
    intercept_from_id = intercept_from_id_;
    slope_from_id = slope_from_id_;
  }
  void insert_step
  (map<double,ThreshIntervals>::iterator it,
   Point point, 
   int high_rank){
    step_ThreshIntervals_map.insert
      (it, 
       pair<double,ThreshIntervals>
       (point.x, ThreshIntervals
        (point.y, high_rank, id_from_rank, rank_from_id)));
  }
  void maybe_add_intersection(double prevStepSize, int high_rank){
    int high_id = (*id_from_rank)[high_rank];
    int low_id = (*id_from_rank)[high_rank-1];
    Point intersectionPoint = intersect
      (intercept_from_id[low_id], intercept_from_id[high_id],
       slope_from_id[low_id], slope_from_id[high_id]);
    // intersection points with infinite values aren't real intersections
    if (intersectionPoint.isFinite() && intersectionPoint.x > prevStepSize) {
      auto it = step_ThreshIntervals_map.lower_bound(intersectionPoint.x);
      if
        (it == step_ThreshIntervals_map.end() || 
         it->first != intersectionPoint.x){
        insert_step(it, intersectionPoint, high_rank);
      }else{
        it->second.add_interval(intersectionPoint.y, high_rank);
      }
    }
  }
};

int lineSearch
(const double *intercept_from_id,
 const double *slope_from_id,
 int lineCount,
 const double *deltaFp,
 const double *deltaFn,
 int maxIterations,
 double *stepSizeVec,
 double *aumVec, 
 double *aumSlopeAfterStepVec,
 double *aucAtStepVec,
 double *aucAfterStepVec,
 int *intersectionCountVec,
 int *intervalCountVec
 ) {
  // map from rank (index when sorted by threshold) to id of line
  // (in FP/FN/etc indices).
  vector<int> id_from_rank(lineCount), rank_from_id(lineCount);
  // map from line number (id) to rank (index in sorted vector) of line.
  vector<double> FP(lineCount), FPlo(lineCount), FPhi(lineCount);
  vector<double> FN(lineCount), FNlo(lineCount), FNhi(lineCount);
  vector<double> M(lineCount);
  for (int i = 0; i < lineCount; i++) {
    id_from_rank[i] = rank_from_id[i] = i;
  }
  Queue queue
    (&id_from_rank, &rank_from_id, intercept_from_id, slope_from_id);
  // start by queueing intersections of every line and the line after it
  for
    (int lineIndexLowBeforeIntersect = 0;
     lineIndexLowBeforeIntersect < lineCount - 1;
     lineIndexLowBeforeIntersect++) {
    int lineIndexHighBeforeIntersect = lineIndexLowBeforeIntersect + 1;
    if
      (intercept_from_id[lineIndexLowBeforeIntersect] >
       intercept_from_id[lineIndexHighBeforeIntersect]) {
      return ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING;
    }
    if
      ((intercept_from_id[lineIndexLowBeforeIntersect] ==
        intercept_from_id[lineIndexHighBeforeIntersect]) &&
       (slope_from_id[lineIndexLowBeforeIntersect] >=
        slope_from_id[lineIndexHighBeforeIntersect])) {
      return ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS;
    }
    queue.maybe_add_intersection(0, lineIndexHighBeforeIntersect);
  }
  double lastStepSize = 0.0;
  double aum = 0.0;
  // build FP, FN, & M, and compute initial aum slope.
  FN[0]=1;
  FN[lineCount - 1] = -deltaFn[lineCount - 1];
  for (int b = lineCount - 2; b >= 1; b--) {
    FN[b] = FN[b + 1] - deltaFn[b];
  }
  FP[0]=0;
  for (int b = 1; b < lineCount; b++) {
    FP[b] = FP[b - 1] + deltaFp[b - 1];
    M[b] = min(FP[b], FN[b]);
    aum += M[b]*(intercept_from_id[b]-intercept_from_id[b-1]);
  }
  double maxFP=0, maxFN=0;
  for(int b=0; b<lineCount; b++){
    maxFP += deltaFp[b];
    maxFN -= deltaFn[b];
  }
  if(maxFP <= 0)return ERROR_LINE_SEARCH_MAX_FP_SHOULD_BE_POSITIVE;
  if(maxFN <= 0)return ERROR_LINE_SEARCH_MAX_FN_SHOULD_BE_POSITIVE;
  // initialize AUC.
  TotalAUC total_auc
    (&FP, maxFP, &FN, maxFN, 
     &M, &id_from_rank, &rank_from_id, slope_from_id, lineCount);
  int last_line = 0;
  double last_thresh = intercept_from_id[0];
  for(int line_i=1; line_i<=lineCount; line_i++){
    double thresh;
    if(line_i==lineCount){
      thresh = INFINITY;
    }else{
      thresh = intercept_from_id[line_i];
    }
    if(last_thresh < thresh){
      total_auc.value += total_auc.get_auc(last_line, line_i);
      last_line = line_i;
      last_thresh = thresh;
    }
  }
  aucAtStepVec[0] = total_auc.value;
  total_auc.zero_update_all();
  // AUM at step size 0
  aumVec[0] = aum;
  aumSlopeAfterStepVec[0] = total_auc.aum_slope;
  stepSizeVec[0] = 0.0;
  aucAfterStepVec[0] = total_auc.value;
  intersectionCountVec[0] = 0;
  intervalCountVec[0] = 0;
  for//iterations/step sizes
    (int iteration = 1; 
     iteration < maxIterations && !queue.step_ThreshIntervals_map.empty(); 
     iteration++){
    auto TI_it = queue.step_ThreshIntervals_map.begin();
    double stepSize = TI_it->first;
    aum += total_auc.aum_slope * (stepSize - lastStepSize);
    stepSizeVec[iteration] = stepSize;
    aumVec[iteration] = aum;
    ThreshIntervals TI = TI_it->second;
    TI.set_intervals_ranks();
    double more_auc_at_step = total_auc.action(&TI, -1.0);
    double auc_after_remove = total_auc.value;
    intersectionCountVec[iteration] = TI.thresh_intervals_map.size();
    intervalCountVec[iteration] = 0;
    for//intersection points. (tie-breaking type 1)
      (auto intervals_it = TI.thresh_intervals_map.begin(); 
       intervals_it != TI.thresh_intervals_map.end();
       intervals_it++){
      double FPhi_tot=0, FPlo_tot=0, FNhi_tot=0, FNlo_tot=0;
      intervalCountVec[iteration] += intervals_it->second.n_intervals;
      int lowest_rank = intervals_it->second.low_rank;
      int highest_rank = intervals_it->second.high_rank;
      for//adjacent lines in an intersection point. (tie-breaking type 2)
        (int low_rank=lowest_rank; 
         low_rank<highest_rank; 
         low_rank++){
        int offset = low_rank-lowest_rank;
        int high_rank = low_rank+1;
        int low_id = id_from_rank[low_rank];
        int top_high_rank = highest_rank-offset;
        int top_high_id = id_from_rank[top_high_rank];
        FPlo_tot += deltaFp[low_id];
        FPlo[high_rank] = FPlo_tot;
        FNlo_tot += deltaFn[low_id];
        FNlo[top_high_rank] = FNlo_tot;
        FPhi_tot += deltaFp[top_high_id];
        FPhi[high_rank] = FPhi_tot;
        FNhi_tot += deltaFn[top_high_id];
        FNhi[top_high_rank] = FNhi_tot;
      }
      for//adjacent lines in an intersection point. (tie-breaking type 2)
        (int low_rank=lowest_rank; 
         low_rank<highest_rank; 
         low_rank++){
        // (∆FP of top line) - (∆FP of bottom line)
        int high_rank = low_rank+1;
        double deltaFpDiff = FPhi[high_rank] - FPlo[high_rank];
        double deltaFnDiff = FNhi[high_rank] - FNlo[high_rank];
        FP[high_rank] += deltaFpDiff;
        FN[high_rank] += deltaFnDiff;
        M[high_rank] = min(FP[high_rank], FN[high_rank]);
      }
      reverse
        (id_from_rank.begin()+lowest_rank, 
         id_from_rank.begin()+highest_rank+1);
      for(int rank = lowest_rank; rank <= highest_rank; rank++){
        int id = id_from_rank[rank];
        rank_from_id[id] = rank;
      }
    }
    total_auc.action(&TI, 1.0);
    aumSlopeAfterStepVec[iteration] = total_auc.aum_slope;
    aucAtStepVec[iteration] = auc_after_remove+more_auc_at_step;
    aucAfterStepVec[iteration] = total_auc.value;
    // queue the next actions/intersections.
    queue.step_ThreshIntervals_map.erase(TI_it);
    int prev_high_rank = 0;
    for//intersection points. (tie-breaking type 1)
      (auto intervals_it = TI.thresh_intervals_map.begin(); 
       intervals_it != TI.thresh_intervals_map.end();
       intervals_it++){
      int highest_rank = intervals_it->second.high_rank;
      int higherRank = highest_rank + 1;
      if (higherRank < lineCount) {
        queue.maybe_add_intersection(stepSize, higherRank);
      }
      int lowest_rank = intervals_it->second.low_rank;
      if (lowest_rank > prev_high_rank) {
        queue.maybe_add_intersection(stepSize, lowest_rank);
      }
      prev_high_rank = highest_rank;
    }
    lastStepSize = stepSize;
  }
  return 0;//SUCCESS
}
