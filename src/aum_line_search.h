#ifndef AUM_LINE_SEARCH_AUMLINESEARCH_H
#define AUM_LINE_SEARCH_AUMLINESEARCH_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
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

/**
 * Line Search
 * @param lines lines array
 * @param lineCount size of the lines array
 * @param deltaFp input FP, assumed length == lineCount
 * @param deltaFn input FN, assumed length == lineCount
 * @param FP destination for FP, assumed length == lineCount
 * @param FN destination for FN, assumed length == lineCount
 * @param M destination for M, assumed length == lineCount
 * @param maxIterations max iterations of the algorithm
 * @param stepSizeVec output vector for the step size of the AUM
 * @param aumVec output vector for the AUM values at the corresponding step sizes
 */
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
);

#define ERROR_LINE_SEARCH_INTERCEPTS_SHOULD_BE_NON_DECREASING 1
#define ERROR_LINE_SEARCH_SLOPES_SHOULD_BE_INCREASING_FOR_EQUAL_INTERCEPTS 2

#endif //AUM_LINE_SEARCH_AUMLINESEARCH_H
