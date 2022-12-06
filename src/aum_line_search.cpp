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

int lineSearch(
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
) {
    // a list of indices of lines
    // this is updated as we move across the X axis, passing intersection points changes the indices
    vector<int> sortedIndices;
    // map from line number to index of line
    vector<int> backwardsIndices;
    sortedIndices.reserve(lineCount);
    backwardsIndices.reserve(lineCount);
    for (int i = 0; i < lineCount; i++) {
        sortedIndices.push_back(i);
        backwardsIndices.push_back(i);
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

    FN[lineCount - 1] = -deltaFn[lineCount - 1];
    for (int b = lineCount - 2; b >= 1; b--) {
        FN[b] = FN[b + 1] - deltaFn[b];
    }

    // build FP, FN, & M, and compute initial aum slope
    for (int b = 1; b < lineCount; b++) {
        double slopeDiff = lines[b].slope - lines[b - 1].slope;
        FP[b] = FP[b - 1] + deltaFp[b - 1];
        M[b] = min(FP[b], FN[b]);
        aumSlope += slopeDiff * M[b];
    }
    double aum = initialAum;

    for (int iterations = 1; iterations < maxIterations && !intersections.empty(); iterations++) {
        auto intersection = *intersections.begin();
        // swap the indices of the lines that make this intersection point
        // the names of these variables look backwards because they get swapped
        int lineIndexLowAfterIntersect = backwardsIndices[intersection.lineLowBeforeIntersect];
        int lineIndexHighAfterIntersect = backwardsIndices[intersection.lineHighBeforeIntersect];
        swap(backwardsIndices[intersection.lineLowBeforeIntersect], backwardsIndices[intersection.lineHighBeforeIntersect]);
        swap(sortedIndices[lineIndexHighAfterIntersect], sortedIndices[lineIndexLowAfterIntersect]);

        // indices of the next lines we want to find intersections for
        int higherLineIndex = lineIndexHighAfterIntersect + 1;
        int lowerLineIndex = lineIndexLowAfterIntersect - 1;

        // (∆FP of top line) - (∆FP of bottom line)
        double deltaFpDiff = deltaFp[intersection.lineHighBeforeIntersect] - deltaFp[intersection.lineLowBeforeIntersect];
        double deltaFnDiff = deltaFn[intersection.lineHighBeforeIntersect] - deltaFn[intersection.lineLowBeforeIntersect];

        // b ∈ {2, . . . , B} is the index of the function which is larger before intersection point
        int b = lineIndexHighAfterIntersect;
        // current step size we're at
        double stepSize = intersection.point.x;

        // update FP & FN
        FP[b] += deltaFpDiff;
        FN[b] += deltaFnDiff;
        double minBeforeIntersection = M[b];
        M[b] = min(FP[b], FN[b]);

        // queue the next intersections in the multiset
        // this creates an intersection between "lineIndexHighAfterIntersect" and "higherLineIndex"
        // "lineIndexHighAfterIntersect" will now be the index of the low line before this new intersection point
        if (higherLineIndex < lineCount) {
            queueIntersection(stepSize, lines, intersections,
                              intersection.lineLowBeforeIntersect, sortedIndices[higherLineIndex]);
        }
        if (lowerLineIndex >= 0) {
            queueIntersection(stepSize, lines, intersections,
                              sortedIndices[lowerLineIndex], intersection.lineHighBeforeIntersect);
        }

        double aumDiff = aumSlope * (stepSize - lastStepSize);
        aum += aumDiff;

        stepSizeVec[iterations] = stepSize;
        aumVec[iterations] = aum;

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
