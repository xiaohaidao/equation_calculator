
#include "Math.h"

#include <cmath>

bool IsEqual(double a, double b) {
    return fabs(a - b) < 1E-12;
}

bool IsGreater(double a, double b) {
    return !IsEqual(a, b) && a > b;
}

bool IsLesser(double a, double b) {
    return !IsEqual(a, b) && a < b;
}

bool IsZero(double x) {
    return IsEqual(x, 0);
}

