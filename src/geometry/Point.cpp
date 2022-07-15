#include "Objects.h"
#include "compare_double.h"
#include "quadratic.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>

Point::Point() {

    x = 0.0;
    y = 0.0;
    z = 0.0;
}

Point::Point(double _x, double _y, double _z) {

    x = _x;
    y = _y;
    z = _z;
}

Point::Point(double _x, double _y) {

    x = _x;
    y = _y;
    z = 0.0;
}
Point::~Point() {}
