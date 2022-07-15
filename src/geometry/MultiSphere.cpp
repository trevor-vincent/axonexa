#include "Objects.h"
#include "compare_double.h"
#include "log.h"
#include "quadratic.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

MultiSphere::MultiSphere(double center_x, double center_y, double center_z,
                         double radius_smallSphere, double radius_bigSphere,
                         double T2_smallSphere, double T2_bigSphere,
                         double D_smallSphere, double D_bigSphere,
                         int region_smallSphere, int region_bigSphere) {

    smallSphere = Sphere(center_x, center_y, center_z, radius_smallSphere,
                         T2_smallSphere, D_smallSphere, region_smallSphere);
    bigSphere = Sphere(center_x, center_y, center_z, radius_bigSphere,
                       T2_bigSphere, D_bigSphere, region_bigSphere);
}

MultiSphere::MultiSphere(Point &r, double radius_smallSphere,
                         double radius_bigSphere, double T2_smallSphere,
                         double T2_bigSphere, double D_smallSphere,
                         double D_bigSphere, int region_smallSphere,
                         int region_bigSphere) {

    smallSphere = Sphere(r.x, r.y, r.z, radius_smallSphere, T2_smallSphere,
                         D_smallSphere, region_smallSphere);
    bigSphere = Sphere(r.x, r.y, r.z, radius_bigSphere, T2_bigSphere,
                       D_bigSphere, region_bigSphere);
}

bool MultiSphere::inside(Vector3 &r) { return bigSphere.inside(r); }

int MultiSphere::getregion(Vector3 &r) {

    int region1 = smallSphere.getregion(r);
    int region2 = bigSphere.getregion(r);

    if (region1 != 0) {

        return region1;

    }

    else if (region2 != 0) {

        return region2;
    }

    return 0;
}

double MultiSphere::getT2(Vector3 &r) {

    double T2_1 = smallSphere.getT2(r);
    double T2_2 = bigSphere.getT2(r);

    if (T2_1 > 0.0) {

        return T2_1;

    }

    else if (T2_2 > 0.0) {

        return T2_2;
    }

    return -1.0;
}

double MultiSphere::getD(Vector3 &r) {

    double D_1 = smallSphere.getD(r);
    double D_2 = bigSphere.getD(r);

    if (D_1 > 0.0) {

        return D_1;

    }

    else if (D_2 > 0.0) {

        return D_2;
    }

    return -1.0;
}

bool MultiSphere::intersect(Line &line, double &v) {

    std::vector<double> pos_v;
    double v1, v2;

    if (smallSphere.intersect(line, v1)) {
        pos_v.push_back(v1);
    }
    if (bigSphere.intersect(line, v2)) {
        pos_v.push_back(v2);
    }

    // FILE_LOG(logDEBUG4) << "*********MULTISPHERE INTERSECTION
    // METHOD*********** " << std::endl; FILE_LOG(logDEBUG4) << "v1 = " << v1 <<
    // std::endl; FILE_LOG(logDEBUG4) << "v2 = " << v2 << std::endl;

    if (pos_v.empty()) {
        return false;
    }

    std::sort(pos_v.begin(), pos_v.end());
    v = pos_v[0];
    return true;
}

Vector3 MultiSphere::getNormal(double x, double y, double z) {
    return smallSphere.getNormal(x, y, z);
}
