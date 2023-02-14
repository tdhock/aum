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

class Actions {
  public:
  void print(){
    for(auto it=ranks.begin(); it != ranks.end(); it++){
      printf("%d,%d\n", it->first, it->second);
    }
  }
  map<int,int> ranks;
  Actions(int high_rank){
    insert_pair(high_rank-1, high_rank);
  }
  void insert_pair(int low, int high){
    ranks.insert(pair<int,int>(low, high));
  }
  void insert_hint_pair(map<int,int>::iterator hint, int low, int high){
    ranks.insert(hint, pair<int,int>(low, high));
  }
  void add_interval(int high_rank){
    int low_rank = high_rank-1;
    auto it_after = ranks.lower_bound(high_rank);
    if(it_after != ranks.begin()){
      auto it_before = it_after;
      it_before--;
      printf("before %d %d\n", it_before->first, it_before->second);
      if(it_before->second == low_rank){
        it_before->second = high_rank;
        return;
      }
    }
    if(it_after != ranks.end()){
      printf("after %d %d\n", it_after->first, it_after->second);
      if(it_after->first == high_rank){
        int new_high = it_after->second;
        auto hint = ranks.erase(it_after);
        insert_hint_pair(hint, low_rank, new_high);
        return;
      }
    }
    insert_hint_pair(it_after, low_rank, high_rank);
  }
};

class TotalAUC {
  public:
  vector<double> *FPR, *FNR, *M;
  vector<int> *id_from_rank;
  const Line *lines;
  double value, aum_slope;
  int lineCount;
  TotalAUC
  (vector<double> *FPR_, 
   vector<double> *FNR_, 
   vector<double> *M_, 
   vector<int> *id_from_rank_,
   const Line *lines_,
   int lineCount_){
    FPR = FPR_;
    FNR = FNR_;
    M = M_;
    id_from_rank = id_from_rank_;
    lines = lines_;
    lineCount = lineCount_;
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
    //printf("%d %f,%f -> %d %f,%f\n", leftRank, fpr(leftRank), tpr(leftRank), rightRank, fpr(rightRank), tpr(rightRank));
    double FPR_diff = fpr(rightRank)-fpr(leftRank);
    double TPR_sum = tpr(rightRank)+tpr(leftRank);
    return FPR_diff*TPR_sum/2;
  }
  void update
  (int first_rightRank, 
   int last_rightRank,
   double sign){
    //printf("update(%d,%d,%f)\n", first_rightRank,  last_rightRank,  sign);
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
      aum_slope += sign*rank_slope*min_diff;
    }
  }
  double action(Actions *act_ptr, double sign){
    double more_auc_at_step = 0;
    for//intersection points. (tie-breaking type 1)
      (auto ranks_it = act_ptr->ranks.begin(); 
       ranks_it != act_ptr->ranks.end();
       ranks_it++){
      int low_rank=ranks_it->first+1;
      int high_rank=ranks_it->second+1;
      //printf("action low=%d hi=%d sign=%f\n", low_rank, high_rank, sign);
      update(low_rank, high_rank, sign);
      more_auc_at_step += get_auc(ranks_it->first, high_rank);
    }
    return more_auc_at_step;
  }
};

class Queue {
  public:
  map<double,Actions> actions;
  vector<int> *id_from_rank;
  const Line *lines;
  int iteration;
  void print(){
    for(auto it=actions.begin(); it != actions.end(); it++){
      printf("step=%f\n", it->first);
      it->second.print();
    }
  }
  Queue(vector<int> *id_from_rank_, const Line *lines_){
    id_from_rank = id_from_rank_;
    lines = lines_;
  }
  void add_intersection(double prevStepSize, int high_rank){
    int high_id = (*id_from_rank)[high_rank];
    int low_id = (*id_from_rank)[high_rank-1];
    auto intersectionPoint = intersect(lines[low_id], lines[high_id]);
    // intersection points with infinite values aren't real intersections
    if (intersectionPoint.isFinite() && intersectionPoint.x > prevStepSize) {
      auto it = actions.lower_bound(intersectionPoint.x);
      //printf("after lower_bound high_rank=%d size=%d step=%f\n", high_rank, actions.size(), intersectionPoint.x);
      if(it == actions.end() || it->first != intersectionPoint.x){
        actions.insert
          (it, pair<double,Actions>(intersectionPoint.x, Actions(high_rank)));
      }else{
        printf("it=%d step.size=%f add_interval(%d) x=%f y=%f\n", iteration, prevStepSize, high_rank, intersectionPoint.x, intersectionPoint.y);
        it->second.add_interval(high_rank);
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
    vector<int> id_from_rank(lineCount);
    // map from line number (id) to rank (index in sorted vector) of line.
    vector<double> FP(lineCount), FPlo(lineCount), FPhi(lineCount);
    vector<double> FN(lineCount), FNlo(lineCount), FNhi(lineCount);
    vector<double> M(lineCount);
    for (int i = 0; i < lineCount; i++) {
        id_from_rank[i] = i;
    }
    Queue queue(&id_from_rank, lines);
    //print_ids(id_from_rank);
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
    TotalAUC total_auc(&FP, &FN, &M, &id_from_rank, lines, lineCount);
    total_auc.value = 0;
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
    total_auc.value = 0;
    total_auc.aum_slope = 0;
    total_auc.update(1, lineCount, 1.0);
    // AUM at step size 0
    aumVec[0] = aum;
    aumSlopeAfterStepVec[0] = total_auc.aum_slope;
    stepSizeVec[0] = 0.0;
    aucAfterStepVec[0] = total_auc.value;
    intersectionCountVec[0] = 0;
    intervalCountVec[0] = 0;
    for//iterations/step sizes
      (int iteration = 1; 
       iteration < maxIterations && !queue.actions.empty(); 
       iteration++){
      queue.iteration=iteration;
      auto act_it = queue.actions.begin();
      double stepSize = act_it->first;
      aum += total_auc.aum_slope * (stepSize - lastStepSize);
      stepSizeVec[iteration] = stepSize;
      aumVec[iteration] = aum;
      Actions *act_ptr = &act_it->second;
      double more_auc_at_step = total_auc.action(act_ptr, -1.0);
      double auc_after_remove = total_auc.value;
      intersectionCountVec[iteration] = act_ptr->ranks.size();
      intervalCountVec[iteration] = 0;
      for//intersection points. (tie-breaking type 1)
        (auto ranks_it = act_ptr->ranks.begin(); 
         ranks_it != act_ptr->ranks.end();
         ranks_it++){
        double FPhi_tot=0, FPlo_tot=0, FNhi_tot=0, FNlo_tot=0;
        intervalCountVec[iteration] += ranks_it->second-ranks_it->first;
        for//adjacent lines in an intersection point. (tie-breaking type 2)
          (int low_rank=ranks_it->first; 
           low_rank<ranks_it->second; 
           low_rank++){
          int offset = low_rank-ranks_it->first;
          int high_rank = low_rank+1;
          int low_id = id_from_rank[low_rank];
          int top_high_rank = ranks_it->second-offset;
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
          (int low_rank=ranks_it->first; 
           low_rank<ranks_it->second; 
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
          (id_from_rank.begin()+ranks_it->first, 
           id_from_rank.begin()+ranks_it->second+1);
      }
      total_auc.action(act_ptr, 1.0);
      aumSlopeAfterStepVec[iteration] = total_auc.aum_slope;
      aucAtStepVec[iteration] = auc_after_remove+more_auc_at_step;
      aucAfterStepVec[iteration] = total_auc.value;
      // queue the next actions/intersections.
      Actions deleted = *act_ptr;
      queue.actions.erase(act_it);
      int prev_high_rank = 0;
      for//intersection points. (tie-breaking type 1)
        (auto ranks_it = deleted.ranks.begin(); 
         ranks_it != deleted.ranks.end();
         ranks_it++){
        int higherRank = ranks_it->second + 1;
        if (higherRank < lineCount) {
            queue.add_intersection(stepSize, higherRank);
        }
        if (ranks_it->first > prev_high_rank) {
            queue.add_intersection(stepSize, ranks_it->first);
        }
        prev_high_rank = ranks_it->second;
      }
      lastStepSize = stepSize;
    }
    return 0;//SUCCESS
}
