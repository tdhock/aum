#include "aum_line_search.h"

using namespace std;

bool Point::operator==(const Point other) const {
    return fabs(x - other.x) < 1e-6 && fabs(y - other.y) < 1e-6;
}

bool Point::isFinite() const {
    return isfinite(x) && isfinite(y);
}

Point intersect(Line a, Line b) {
    if (a.slope == b.slope) return Point{INFINITY, INFINITY};
    double x = (b.intercept - a.intercept) / (a.slope - b.slope);
    double y = a.intercept + a.slope * x;
    return Point{x, y};
}

bool operator<(IntersectionData a, IntersectionData b) {
    return a.point.x < b.point.x;
}

bool operator==(IntersectionData a, IntersectionData b) {
    return a.point == b.point && (
            (a.lineHighBeforeIntersect == b.lineHighBeforeIntersect &&
             a.lineLowBeforeIntersect == b.lineLowBeforeIntersect) ||
            (a.lineHighBeforeIntersect == b.lineLowBeforeIntersect &&
             a.lineLowBeforeIntersect == b.lineHighBeforeIntersect)
    );
}

void queueIntersection(
        double currentStepSize,
        const Line *lines,
        multiset<IntersectionData> &intersections,
        int lowLine,
        int highLine
) {
    auto intersectionPoint = intersect(
            lines[lowLine],
            lines[highLine]
    );

    // intersection points with infinite values aren't real intersections
    if (intersectionPoint.isFinite() && intersectionPoint.x > currentStepSize) {
        auto intersection = IntersectionData{
                intersectionPoint, lowLine, highLine
        };

        intersections.insert(intersection);
    }
}

class TotalAUC {
  public:
  vector<double> *FPR, *FNR;
  double value;
  int lineCount;
  TotalAUC(vector<double> *FPR_, vector<double> *FNR_, int lineCount_){
    FPR = FPR_;
    FNR = FNR_;
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
  }
};

int lineSearch(
        const Line *lines,
        int lineCount,
        const double *deltaFp,
        const double *deltaFn,
        double initialAum,
        int maxIterations,
        double *stepSizeVec,
        double *aumVec, 
        double *aucAtStepVec,
        double *aucAfterStepVec
) {
    // a list of indices of lines map from rank (index when sorted by
    // threshold) to id of line (in FP/FN/etc indices).
    vector<int> id_from_rank(lineCount);
    // map from line number (id) to rank (index in sorted vector) of line.
    vector<int> rank_from_id(lineCount);
    vector<double> FP(lineCount);
    vector<double> FN(lineCount);
    vector<double> M(lineCount);
    for (int i = 0; i < lineCount; i++) {
        id_from_rank[i] = i;
        rank_from_id[i] = i;
    }
    multiset<IntersectionData> intersections;
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
	queueIntersection
	  (0, lines, intersections,
	   lineIndexLowBeforeIntersect,
	   lineIndexHighBeforeIntersect);
    }

    // AUM at step size 0
    double intercept = initialAum;
    aumVec[0] = intercept;
    stepSizeVec[0] = 0.0;
    double aumSlope = 0.0;
    double lastStepSize = 0.0;

    // build FP, FN, & M, and compute initial aum slope.
    FP[0]=0;
    FN[0]=1;
    FN[lineCount - 1] = -deltaFn[lineCount - 1];
    for (int b = lineCount - 2; b >= 1; b--) {
        FN[b] = FN[b + 1] - deltaFn[b];
    }
    for (int b = 1; b < lineCount; b++) {
        double slopeDiff = lines[b].slope - lines[b - 1].slope;
        FP[b] = FP[b - 1] + deltaFp[b - 1];
        M[b] = min(FP[b], FN[b]);
        aumSlope += slopeDiff * M[b];
    }
    // initialize AUC.
    TotalAUC total_auc(&FP, &FN, lineCount);
    total_auc.value = 0;
    total_auc.update(1, lineCount, 1.0);
    aucAtStepVec[0] = aucAfterStepVec[0] = total_auc.value;
    double aum = initialAum;
    for (int iteration = 1; iteration < maxIterations && !intersections.empty(); iteration++) {
        auto intersection = *intersections.begin();
        // (∆FP of top line) - (∆FP of bottom line)
        double deltaFpDiff = 
          deltaFp[intersection.lineHighBeforeIntersect] - 
          deltaFp[intersection.lineLowBeforeIntersect];
        double deltaFnDiff = 
          deltaFn[intersection.lineHighBeforeIntersect] - 
          deltaFn[intersection.lineLowBeforeIntersect];
        // b ∈ {2, . . . , B} is the rank of the function which is
        // larger before intersection point
        int b = rank_from_id[intersection.lineHighBeforeIntersect];
        // current step size we're at
        double stepSize = intersection.point.x;
        total_auc.update(b, b+1, -1.0);
        aucAtStepVec[iteration] = 
          total_auc.value+total_auc.get_auc(b-1, b+1);
        // update FP & FN
        FP[b] += deltaFpDiff;
        FN[b] += deltaFnDiff;
        double minBeforeIntersection = M[b];
        M[b] = min(FP[b], FN[b]);
        total_auc.update(b, b+1, 1.0);
        aucAfterStepVec[iteration] = total_auc.value;
        swap
          (rank_from_id[intersection.lineLowBeforeIntersect], 
           rank_from_id[intersection.lineHighBeforeIntersect]);
        swap
          (id_from_rank[b], 
           id_from_rank[b-1]);
        // queue the next intersections in the multiset.
        int higherRank = b + 1;
        if (higherRank < lineCount) {
            queueIntersection
              (stepSize, lines, intersections,
               intersection.lineLowBeforeIntersect, 
               id_from_rank[higherRank]);
        }
        int lowerRank = b - 2;
        if (lowerRank >= 0) {
            queueIntersection
              (stepSize, lines, intersections,
               id_from_rank[lowerRank], 
               intersection.lineHighBeforeIntersect);
        }

        double aumDiff = aumSlope * (stepSize - lastStepSize);
        aum += aumDiff;

        stepSizeVec[iteration] = stepSize;
        aumVec[iteration] = aum;

        // update aum slope
        double slopeDiff = lines[intersection.lineHighBeforeIntersect].slope - lines[intersection.lineLowBeforeIntersect].slope;
        // this is the D^(k+1) update rule in the paper,
        // it updates the AUM slope for the next iteration
        double mAfter = b + 1 < lineCount ? M[b + 1] : 0;
        aumSlope += (slopeDiff) * (mAfter + M[b - 1] - M[b] - minBeforeIntersection);

        intersections.erase(intersection);
        lastStepSize = stepSize;
    }
    return 0;//SUCCESS
}
