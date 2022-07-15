
//#include <cmath>
//#include <iostream>
//#include <limits>
//#include <vector>
//#include <algorithm>
//#include <utility>
//#include "Objects.h"
//#include "particle.h"
//#include "ran.h"

#define PI 3.1415926535897932384626433832795

// using namespace std;

#ifndef LATTICE_H_ // only allows one class declaration memory allotment to be
                   // made
#define LATTICE_H_

template <class OBJECT>
bool sortPair(const pair<double, OBJECT> &a, const pair<double, OBJECT> &b) {
    return a.first < b.first;
}

template <class OBJECT> class Lattice {

  private:
    double D_extra;
    double T2_extra;
    double permeability;
    vector<OBJECT> basis;

    // lattice vectors
    double a, b, c;
    Vector3 ahat, bhat, chat;

  public:
    // Lattice(double _D_extra, double _T2_extra, double _permeability, int
    // _number_of_elements){

    // D_extra = _D_extra;
    // T2_extra = _T2_extra;
    // permeability = _permeability;
    // basis.resize(_number_of_elements);

    // }

    Lattice(double _D_extra, double _T2_extra, double _permeability) {

        D_extra = _D_extra;
        T2_extra = _T2_extra;
        permeability = _permeability;
    }

    ~Lattice() {}

    double get_permeability() { return permeability; }

    void initializeUniformly(Ran &generator, vector<Particle> &ensemble) {

        FILE_LOG(logINFO) << " PARTICLES UNIFORMLY DISTRIBUTED ABOUT SUBSTRATE"
                          << endl;
        Point scoords;
        double a0, b0, c0;

        for (int i = 0; i < ensemble.size(); i++) {
            a0 = generator.doub() * a;
            b0 = generator.doub() * b;
            c0 = generator.doub() * c;
            scoords = getsimulcoords(a0, b0, c0);
            ensemble[i] = Particle(scoords.x, scoords.y, scoords.z);
        }
    }

    // initialize particles to basis region
    // note that this method assumes intracellular region is defined as: region
    // != 0
    void initializeInBasis(Ran &generator, vector<Particle> &ensemble) {

        Point scoords;
        double a0, b0, c0;
        int i = 0;
        FILE_LOG(logINFO) << " PARTICLES INITIALIZED IN BASIS " << endl;

        while (i < ensemble.size()) {
            a0 = generator.doub() * a;
            b0 = generator.doub() * b;
            c0 = generator.doub() * c;
            scoords = getsimulcoords(a0, b0, c0);
            Vector3 svec(scoords.x, scoords.y, scoords.z);
            // cout << i << " " << a0 << " " << b0 << " " << c0 << " " << svec
            // << "  " << inregion(svec) << endl;
            if (inregion(svec) != 0) {
                ensemble[i] = Particle(scoords.x, scoords.y, scoords.z);
            } else {
                i--;
            }
            i++;
        }
    }

    void initializeInRegions(Ran &generator, vector<Particle> &ensemble,
                             vector<int> &regions) {

        Point scoords;
        double a0, b0, c0;
        int i = 0;
        FILE_LOG(logINFO) << " PARTICLES INITIALIZED IN SPECIFIED REGIONS "
                          << endl;

        while (i < ensemble.size()) {
            a0 = generator.doub() * a;
            b0 = generator.doub() * b;
            c0 = generator.doub() * c;
            scoords = getsimulcoords(a0, b0, c0);
            Vector3 svec(scoords.x, scoords.y, scoords.z);
            // cout << i << " " << a0 << " " << b0 << " " << c0 << " " << svec
            // << "  " << inregion(svec) << endl;

            bool inside = false;
            for (int j = 0; j < regions.size(); j++) {
                if (inregion(svec) == regions[j]) {
                    ensemble[i] = Particle(scoords.x, scoords.y, scoords.z);
                    inside = true;
                    break;
                }
            }
            if (inside == false) {
                i--;
            }
            i++;
        }
    }

    Point getsimulcoords(double a0, double b0, double c0) {

        Vector3 xhat(1.0, 0.0, 0.0);
        Vector3 yhat(0.0, 1.0, 0.0);
        Vector3 zhat(0.0, 0.0, 1.0);

        double x = (ahat * a0 + bhat * b0 + chat * c0) * xhat;
        double y = (ahat * a0 + bhat * b0 + chat * c0) * yhat;
        double z = (ahat * a0 + bhat * b0 + chat * c0) * zhat;

        return Point(x, y, z);
    }

    Point getlatticecoords(double x, double y, double z) {

        Vector3 xhat(1.0, 0.0, 0.0);
        Vector3 yhat(0.0, 1.0, 0.0);
        Vector3 zhat(0.0, 0.0, 1.0);

        double a0 = (xhat * x + yhat * y + zhat * z) * ahat;
        double b0 = (xhat * x + yhat * y + zhat * z) * bhat;
        double c0 = (xhat * x + yhat * y + zhat * z) * chat;

        return Point(a0, b0, c0);
    }

    void correctBoundary(Vector3 &r) {

        if (inLatticeCell(r.x, r.y, r.z) == false) {

            Point lcoords = getlatticecoords(r.x, r.y, r.z);

            if (lcoords.x > a) {
                lcoords.x = fmod(lcoords.x, a);
            }
            if (lcoords.x < 0.0) {
                lcoords.x = fmod(lcoords.x, a) + a;
            }

            if (lcoords.y > b) {
                lcoords.y = fmod(lcoords.y, b);
            }
            if (lcoords.y < 0.0) {
                lcoords.y = fmod(lcoords.y, b) + b;
            }

            if (lcoords.z > c) {
                lcoords.z = fmod(lcoords.z, c);
            }
            if (lcoords.z < 0.0) {
                lcoords.z = fmod(lcoords.z, c) + c;
            }

            lcoords = getsimulcoords(lcoords.x, lcoords.y, lcoords.z);

            r.x = lcoords.x;
            r.y = lcoords.y;
            r.z = lcoords.z;
        }
    }

    void lmod(Point &r) {

        if (!inLatticeCell(r.x, r.y, r.z)) {

            Point lcoords = getlatticecoords(r.x, r.y, r.z);

            if (lcoords.x > a) {
                lcoords.x = fmod(lcoords.x, a);
            }
            if (lcoords.x < 0.0) {
                lcoords.x = fmod(lcoords.x, a) + a;
            }

            if (lcoords.y > b) {
                lcoords.y = fmod(lcoords.y, b);
            }
            if (lcoords.y < 0.0) {
                lcoords.y = fmod(lcoords.y, b) + b;
            }

            if (lcoords.z > c) {
                lcoords.z = fmod(lcoords.z, c);
            }
            if (lcoords.z < 0.0) {
                lcoords.z = fmod(lcoords.z, c) + c;
            }

            r = getsimulcoords(lcoords.x, lcoords.y, lcoords.z);
        }
    }

    double getT2(Vector3 &x) {

        double T2;

        for (int i = 0; i < basis.size(); i++) {

            T2 = basis[i].getT2(x);
            if (T2 > 0) {
                return T2;
            }
        }

        return T2_extra;
    }

    double getD(Vector3 &x) {

        double D;

        for (int i = 0; i < basis.size(); i++) {

            D = basis[i].getD(x);
            if (D > 0) {
                return D;
            }
        }

        return D_extra;
    }

    int inregion(Vector3 &x) {

        int region;

        for (int i = 0; i < basis.size(); i++) {
            region = basis[i].getregion(x);
            if (region != 0) {
                return region;
            }
        }

        return 0;
    }

    bool inLatticeCell(double x, double y, double z) {

        Point lcoords = getlatticecoords(x, y, z);
        if (lcoords.x > a || lcoords.x < 0.0 || lcoords.y > b ||
            lcoords.y < 0.0 || lcoords.z > c || lcoords.z < 0.0) {
            return false;
        }
        return true;
    }

    bool intersection(Vector3 &final_position, Vector3 &initial_position,
                      Vector3 &normal, double &v) {

        vector<pair<double, Vector3>> pos_v;
        FILE_LOG(logDEBUG4) << endl
                            << " **************** LATTICE INTERSECTION METHOD "
                               "***************** "
                            << endl;
        FILE_LOG(logDEBUG4) << final_position << endl;
        FILE_LOG(logDEBUG4) << initial_position << endl;
        FILE_LOG(logDEBUG4) << "basis size = " << basis.size() << endl;

        Point initial(initial_position.x, initial_position.y,
                      initial_position.z);
        Point final(final_position.x, final_position.y, final_position.z);
        Vector3 dr = final_position - initial_position;

        lmod(initial);
        lmod(final);
        Line line1(final.x - dr.x, final.y - dr.y, final.z - dr.z, final.x,
                   final.y, final.z);
        Line line2(initial.x, initial.y, initial.z, initial.x + dr.x,
                   initial.y + dr.y, initial.z + dr.z);

        for (int i = 0; i < basis.size(); i++) {
            // FILE_LOG(logDEBUG4) << i << endl;
            if (basis[i].intersect(line1, v)) {
                line1 *= v;
                Vector3 norm = basis[i].getNormal(line1.xf, line1.yf, line1.zf);
                pos_v.push_back(std::make_pair(v, norm));
            }
            if (basis[i].intersect(line2, v)) {
                line2 *= v;
                Vector3 norm = basis[i].getNormal(line2.xf, line2.yf, line2.zf);
                pos_v.push_back(std::make_pair(v, norm));
            }
        }

        if (pos_v.empty()) {
            return false;
        }

        sort(pos_v.begin(), pos_v.end(), sortPair<Vector3>);
        v = pos_v[0].first;
        normal = pos_v[0].second;

        FILE_LOG(logDEBUG4)
            << endl
            << " ********************************************* " << endl;

        return true;
    }

    void addBasis(OBJECT o) { basis.push_back(o); }

    void setLatticeVectors(double _a, double _b, double _c, Vector3 _ahat,
                           Vector3 _bhat, Vector3 _chat) {

        a = _a;
        b = _b;
        c = _c;
        ahat = _ahat;
        bhat = _bhat;
        chat = _chat;
    }

    void setLatticeVectors(Vector3 _a, Vector3 _b, Vector3 _c) {
        a = _a.magnitude();
        b = _b.magnitude();
        c = _c.magnitude();
        _a.normalize();
        _b.normalize();
        _c.normalize();
        ahat = _a;
        bhat = _b;
        chat = _c;
    }

    double calcf(int points, int seed) {

        Ran randomnumber(seed);
        Point rand_simul;
        double temp[3];
        int points_inside = 0;

        for (int i = 0; i < points; i++) {

            rand_simul =
                getsimulcoords(a * randomnumber.doub(), b * randomnumber.doub(),
                               c * randomnumber.doub());
            temp[0] = rand_simul.x;
            temp[1] = rand_simul.y;
            temp[2] = rand_simul.z;

            for (int j = 0; j < basis.size(); j++) {

                if (basis[j].inside(temp)) {
                    points_inside++;
                }
            }
        }

        return ((double)points_inside) / ((double)points);
    }
};

#endif
