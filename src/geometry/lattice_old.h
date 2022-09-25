
#include "Objects.h"
#include "compare_double.h"
#include "quadratic.h"
#include "ran.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#define PI 3.1415926535897932384626433832795

using namespace std;

#ifndef LATTICE_H_ // only allows one class declaration memory allotment to be
                   // made
#define LATTICE_H_

/*

The main parent class Lattice. It serves as only a label so that the user can
choose the specific type of lattice at runtime and allow me to still make
generic routines for only Lattice but still use them for the specific lattice
types.

*/

class SzaferBoxLattice : public Lattice<SquareCylinder_XY> {

  public:
    SzaferBoxLattice(double, double, double, double, double, double, double);
    ~SzaferBoxLattice();
};

template <class OBJECT> class Lattice {

  protected:
    double D_extra;
    double T2_extra;
    double permeability;

    vector<OBJECT> basis;

    // lattice vectors
    double a, b, c;
    Vector3 ahat, bhat, chat;

  public:
    Lattice(double, double, double, int);
    ~Lattice();

    double get_permeability();
    void initialize(Ran, double[][3], int);
    Point getsimulcoords(double, double, double);
    Point getlatticecoords(double, double, double);
    void correctBoundary(double[]);
    void lmod(Line &line);
    double getT2(double[]);
    double getD(double[]);
    int inregion(double[]);
    bool inLatticeCell(double, double, double);
    double intersection(double[], double[], double, double[]);
    void addBasis(OBJECT);
    void setLatticeVectors(double, double, Vector3, Vector3);
};

#endif
