#ifndef AUM_LINE_SEARCH_AUMLINESEARCH_H
#define AUM_LINE_SEARCH_AUMLINESEARCH_H
#include <algorithm>//reverse
#include <vector>
#include <map>
#include <cmath>//isfinite

struct Point {
    double x; // step size
    double y;
    bool operator==(Point other) const;
    bool isFinite() const;
};

int lineSearch(
        const double *intercept,
        const double *slope,
        const int lineCount,
        const double *deltaFp,
        const double *deltaFn,
        const int maxIterations,
        //inputs above, size of outputs below = maxIterations.
        double *stepSizeVec,
        double *aumVec,
        double *aumSlopeAfterStepVec,
        double *aucAtStepVec,
        double *aucAfterStepVec,
        int *intersectionCountVec,
        int *intervalCountVec,
	int *qSizeVec
);

#define ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING 1
#define ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS 2
#define ERROR_LINE_SEARCH_MAX_FP_SHOULD_BE_POSITIVE 3
#define ERROR_LINE_SEARCH_MAX_FN_SHOULD_BE_POSITIVE 3

#endif //AUM_LINE_SEARCH_AUMLINESEARCH_H
