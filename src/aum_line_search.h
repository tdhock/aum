#ifndef AUM_LINE_SEARCH_AUMLINESEARCH_H
#define AUM_LINE_SEARCH_AUMLINESEARCH_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cmath>
#include <unordered_set>

struct Line {
    double slope;
    double intercept;
};

struct Point {
    double x; // step size
    double y;

    bool operator==(Point other) const;

    bool isFinite() const;
};

Point intersect(Line a, Line b);

struct IntersectionData {
    Point point;
    // "low" line is the line below the other before the intersection point
    int lineLowBeforeIntersect;
    // "high" line is the line above the other after the intersection point
    int lineHighBeforeIntersect;
};

namespace std {
    // hash implementation for Point used in the multiset below
    template<>
    struct hash<Point> {
        size_t operator()(const Point &point) const noexcept {
            return ((uint64_t) point.x) << 32 | ((uint64_t) point.y);
        }
    };
}

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
void lineSearch(
        const Line *lines,
        int lineCount,
        const double *deltaFp,
        const double *deltaFn,
        double initialAum,
        int maxIterations,
        double *FP,
        double *FN,
        double *M,
        double *stepSizeVec,
        double *aumVec
);

#endif //AUM_LINE_SEARCH_AUMLINESEARCH_H
