#include "aum_line_search.h"
using namespace std;

bool Point::isFinite() const {
    return isfinite(x) && isfinite(y);
}

Point intersect(Line a, Line b) {
    if (a.slope == b.slope) return Point{INFINITY, INFINITY};
    double x = (b.intercept - a.intercept) / (a.slope - b.slope);
    double y = a.intercept + a.slope * x;
    return Point{x, y};
}

class Intervals {
  public:
  int high_id;
  int n_intervals;
  vector<int> *id_from_rank, *rank_from_id;
  int get_high_rank(){
    return (*rank_from_id)[high_id];
  }
  int get_low_rank(){
    return get_high_rank()-n_intervals;
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
       (thresh, Intervals{high_id,1,id_from_rank,rank_from_id}));
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
          printf("WARNING: in add_interval, old_low_rank=%d old_high_rank=%d new_high_rank=%d which should never happen!\n", old_low_rank, old_high_rank, new_high_rank);
        }
        it_after_or_same->second.n_intervals++;
        return;
      }
    }
    thresh_intervals_map.insert
      (it_after_or_same, 
       pair<double,Intervals>
       (thresh, Intervals{new_high_id,1,id_from_rank,rank_from_id}));
  }
};

class TotalAUC {
  public:
  vector<double> *FPR, *FNR, *M;
  vector<int> *id_from_rank, *rank_from_id;
  const Line *lines;
  double value, aum_slope;
  int lineCount;
  TotalAUC
  (vector<double> *FPR_, 
   vector<double> *FNR_, 
   vector<double> *M_, 
   vector<int> *id_from_rank_,
   vector<int> *rank_from_id_,
   const Line *lines_,
   int lineCount_){
    FPR = FPR_;
    FNR = FNR_;
    M = M_;
    id_from_rank = id_from_rank_;
    rank_from_id = rank_from_id_;
    lines = lines_;
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
    return 1-(*FNR)[rank];
  }
  double fpr(int rank){
    if(rank==lineCount)return 1;
    return (*FPR)[rank];
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
      double rank_slope = lines[id].slope;
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
      int high_rank=(*rank_from_id)[TI_it->second.high_id]+1;
      int low_rank=high_rank-TI_it->second.n_intervals;
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
  const Line *lines;
  int iteration;
  Queue
  (vector<int> *id_from_rank_, 
   vector<int> *rank_from_id_,
   const Line *lines_){
    id_from_rank = id_from_rank_;
    rank_from_id = rank_from_id_;
    lines = lines_;
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
  void add_intersection(double prevStepSize, int high_rank){
    int high_id = (*id_from_rank)[high_rank];
    int low_id = (*id_from_rank)[high_rank-1];
    Point intersectionPoint = intersect(lines[low_id], lines[high_id]);
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

int lineSearch(
        const Line *lines,
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
    Queue queue(&id_from_rank, &rank_from_id, lines);
    // start by queueing intersections of every line and the line after it
    for (int lineIndexLowBeforeIntersect = 0;
	 lineIndexLowBeforeIntersect < lineCount - 1;
	 lineIndexLowBeforeIntersect++) {
        int lineIndexHighBeforeIntersect = lineIndexLowBeforeIntersect + 1;
	if (lines[lineIndexLowBeforeIntersect].intercept >
	    lines[lineIndexHighBeforeIntersect].intercept) {
	  return ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING;
	}
	if ((lines[lineIndexLowBeforeIntersect].intercept ==
	     lines[lineIndexHighBeforeIntersect].intercept) &&
	    (lines[lineIndexLowBeforeIntersect].slope >=
	     lines[lineIndexHighBeforeIntersect].slope)) {
	  return ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS;
	}
	queue.add_intersection(0, lineIndexHighBeforeIntersect);
    }
    double lastStepSize = 0.0;
    double aum = 0.0;
    // build FP, FN, & M, and compute initial aum slope.
    FP[0]=0;
    FN[0]=1;
    FN[lineCount - 1] = -deltaFn[lineCount - 1];
    for (int b = lineCount - 2; b >= 1; b--) {
        FN[b] = FN[b + 1] - deltaFn[b];
    }
    for (int b = 1; b < lineCount; b++) {
        FP[b] = FP[b - 1] + deltaFp[b - 1];
        M[b] = min(FP[b], FN[b]);
        aum += M[b]*(lines[b].intercept-lines[b-1].intercept);
    }
    // initialize AUC.
    TotalAUC total_auc
      (&FP, &FN, &M, &id_from_rank, &rank_from_id, lines, lineCount);
    int last_line = 0;
    double last_thresh = lines[0].intercept;
    for(int line_i=1; line_i<=lineCount; line_i++){
      double thresh;
      if(line_i==lineCount){
        thresh = INFINITY;
      }else{
        thresh = lines[line_i].intercept;
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
      queue.iteration=iteration;
      auto TI_it = queue.step_ThreshIntervals_map.begin();
      double stepSize = TI_it->first;
      aum += total_auc.aum_slope * (stepSize - lastStepSize);
      stepSizeVec[iteration] = stepSize;
      aumVec[iteration] = aum;
      ThreshIntervals *TI_ptr = &TI_it->second;
      double more_auc_at_step = total_auc.action(TI_ptr, -1.0);
      double auc_after_remove = total_auc.value;
      intersectionCountVec[iteration] = TI_ptr->thresh_intervals_map.size();
      intervalCountVec[iteration] = 0;
      for//intersection points. (tie-breaking type 1)
        (auto intervals_it = TI_ptr->thresh_intervals_map.begin(); 
         intervals_it != TI_ptr->thresh_intervals_map.end();
         intervals_it++){
        double FPhi_tot=0, FPlo_tot=0, FNhi_tot=0, FNlo_tot=0;
        intervalCountVec[iteration] += intervals_it->second.n_intervals;
        int lowest_rank = intervals_it->second.get_low_rank();
        int highest_rank = intervals_it->second.get_high_rank();
        for//adjacent lines in an intersection point. (tie-breaking type 2)
          (int low_rank=lowest_rank; 
           low_rank<highest_rank; 
           low_rank++){
          int offset = low_rank-intervals_it->first;
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
      total_auc.action(TI_ptr, 1.0);
      aumSlopeAfterStepVec[iteration] = total_auc.aum_slope;
      aucAtStepVec[iteration] = auc_after_remove+more_auc_at_step;
      aucAfterStepVec[iteration] = total_auc.value;
      // queue the next actions/intersections.
      ThreshIntervals deleted = *TI_ptr;
      queue.step_ThreshIntervals_map.erase(TI_it);
      int prev_high_rank = 0;
      for//intersection points. (tie-breaking type 1)
        (auto intervals_it = deleted.thresh_intervals_map.begin(); 
         intervals_it != deleted.thresh_intervals_map.end();
         intervals_it++){
        int higherRank = intervals_it->second + 1;
        if (higherRank < lineCount) {
            queue.add_intersection(stepSize, higherRank);
        }
        if (intervals_it->first > prev_high_rank) {
            queue.add_intersection(stepSize, intervals_it->first);
        }
        prev_high_rank = intervals_it->second;
      }
      lastStepSize = stepSize;
    }
    return 0;//SUCCESS
}
