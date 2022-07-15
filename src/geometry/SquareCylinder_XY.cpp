#include "Objects.h"
#include "compare_double.h"
#include "quadratic.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>

SquareCylinder_XY::SquareCylinder_XY(double _center_x, double _center_y,
                                     double _width, int _region, double _T2,
                                     double _D) {

    center_x = _center_x;
    center_y = _center_y;
    width = _width;
    region = _region;
    T2 = _T2;
    D = _D;

    x_corner1 = center_x - width / 2.0;
    y_corner1 = center_y - width / 2.0;
    x_corner2 = center_x + width / 2.0;
    y_corner2 = center_y + width / 2.0;
}

bool SquareCylinder_XY::inside(Vector3 &r) {

    if ((r.x > (center_x - width / 2.0)) && (r.x < (center_x + width / 2.0)) &&
        (r.y > center_y - width / 2.0) && (r.y < center_y + width / 2.0)) {
        return true;
    }

    return false;
}

Vector3 SquareCylinder_XY::getNormal(double x, double y, double z) {

    double n_x = 0.0;
    double n_y = 0.0;

    // std::cout << x << " " << y << " " << z << std::endl;

    if (doub_equal(y, center_y + width / 2.0) == true) {
        n_y = 1.0;
    }
    if (doub_equal(x, center_x + width / 2.0) == true) {
        n_x = 1.0;
    }
    if (doub_equal(y, center_y - width / 2.0) == true) {
        n_y = 1.0;
    }
    if (doub_equal(x, center_x - width / 2.0) == true) {
        n_x = 1.0;
    }

    // std::cout << " center_y + width/2.0 = " << center_y + width/2.0 <<
    // std::endl; std::cout << " center_y - width/2.0 = " << center_y -
    // width/2.0 << std::endl; std::cout << " center_x + width/2.0 = " <<
    // center_x + width/2.0 << std::endl; std::cout << " center_x - width/2.0 =
    // " << center_x - width/2.0 << std::endl;

    // std::cout << " normal = " << n_x << " " << n_y << std::endl;

    return Vector3(n_x, n_y, 0.0);
}

int SquareCylinder_XY::getregion(Vector3 &r) {

    if (inside(r)) {
        return region;
    }

    return 0;
}

double SquareCylinder_XY::getT2(Vector3 &r) {

    if (inside(r)) {
        return T2;
    }

    return -1.0;
}

double SquareCylinder_XY::getD(Vector3 &r) {

    if (inside(r)) {
        return D;
    }

    return -1.0;
}

bool SquareCylinder_XY::intersect(Line &line, double &v) {

    std::vector<double> pos_v;

    double v1, v2, v3, v4, dx, dy, dz;
    double xi = line.xi;
    double xf = line.xf;
    double yi = line.yi;
    double yf = line.yf;
    double zi = line.zi;
    double zf = line.zf;

    dx = xf - xi;
    dy = yf - yi;
    dz = zf - zi;

    double step_mag = sqrt(dx * dx + dy * dy + dz * dz);

    v1 = (y_corner1 - yi) / (yf - yi); // v for face 1
    v2 = (x_corner1 - xi) / (xf - xi); // v for face 2
    v3 = (y_corner2 - yi) / (yf - yi); // v for face 3
    v4 = (x_corner2 - xi) / (xf - xi); // v for face 4

    if (v1 > 0.0 && v1 < 1.0 && v1 * dx + xi < x_corner2 &&
        v1 * dx + xi > x_corner1 && doub_isnan(v1) == false &&
        doub_equal(v1 * step_mag, 0.0) == false) {
        pos_v.push_back(v1);
    }
    if (v2 > 0.0 && v2 < 1.0 && v2 * dy + yi < y_corner2 &&
        v2 * dy + yi > y_corner1 && doub_isnan(v2) == false &&
        doub_equal(v2 * step_mag, 0.0) == false) {
        pos_v.push_back(v2);
    }
    if (v3 > 0.0 && v3 < 1.0 && v3 * dx + xi < x_corner2 &&
        v3 * dx + xi > x_corner1 && doub_isnan(v3) == false &&
        doub_equal(v3 * step_mag, 0.0) == false) {
        pos_v.push_back(v3);
    }
    if (v4 > 0.0 && v4 < 1.0 && v4 * dy + yi < y_corner2 &&
        v4 * dy + yi > y_corner1 && doub_isnan(v4) == false &&
        doub_equal(v4 * step_mag, 0.0) == false) {
        pos_v.push_back(v4);
    }
    if (pos_v.empty()) {
        return false;
    }

    std::sort(pos_v.begin(), pos_v.end());
    v = pos_v[0];
    return true;
}
