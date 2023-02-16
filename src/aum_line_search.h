#ifndef AUM_LINE_SEARCH_AUMLINESEARCH_H
#define AUM_LINE_SEARCH_AUMLINESEARCH_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>

class Line {
  public:
  double slope;
  double intercept;
  double thresh(double step) const;
};

struct Point {
    double x; // step size
    double y;

    bool operator==(Point other) const;

    bool isFinite() const;
};

Point intersect(Line a, Line b);

int lineSearch(
        const Line *lines,
        int lineCount,
        const double *deltaFp,
        const double *deltaFn,
        int maxIterations,
        //inputs above, size of outputs below = maxIterations.
        double *stepSizeVec,
        double *aumVec,
        double *aumSlopeAfterStepVec,
        double *aucAtStepVec,
        double *aucAfterStepVec,
        int *intersectionCountVec,
        int *intervalCountVec
);

#define ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING 1
#define ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS 2

#endif //AUM_LINE_SEARCH_AUMLINESEARCH_H
