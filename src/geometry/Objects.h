

#ifndef OBJECTS_H_ // only allows one class declaration memory allotment to be
                   // made
#define OBJECTS_H_

//#include "compare_double.h"

class Line {

  public:
    double xi;
    double yi;

    double xf;
    double yf;

    double zi;
    double zf;

    Line(double, double, double, double, double, double);
    Line(double, double, double, double);
    double magnitude();
    void operator*=(const double);
    Line operator*(const double) const;
    ~Line();
};

class Point {

  public:
    double x;
    double y;
    double z;

    Point();
    Point(double, double, double);
    Point(double, double);
    ~Point();
};

/* Holds a vector in three dimensions. Four data members are allocated to ensure
alignment in an array. */

class Vector3 {

  public:
    /* Holds the value along the x axis. */
    double x;

    /* Holds the value along the y axis */
    double y;

    /* Holds the value along the z axis */
    double z;

  private:
    /* I'v added an extra piece of data in the vector structure
            called pad. This isn't part of the mathematics of vectors, and is
       purely there for performance. On many machines, four floating-point
       values sit more cleanly in memory than three (memory is optimized for
       sets of four words), so noticeable speed-ups can be achieved by adding
       this padding.
    */
    double pad;

  public:
    /* The default constructor creates a zero vector */
    Vector3() : x(0), y(0), z(0) {}

    /* The explicit constructor creates a vector with the
    given components		*/

    Vector3(const double x, const double y, const double z)
        : x(x), y(y), z(z) {}

    /* Flips all the components of the vector */

    void invert() {

        x = -x;
        y = -y;
        z = -z;
    }

    /* Gets the magnitude of this vector. */
    double magnitude() const { return sqrt(x * x + y * y + z * z); }

    /* Gets the squared magnitude of this vector. In some
    cases we do not need the exact magnitude; for example, when
    we need to compare two magnitudes to see which is greater, it
    is faster the compare the squares of the magnitudes. This is
    b/c we do not need to call the sqrt function.*/
    double squareMagnitude() const { return x * x + y * y + z * z; }

    /* Turns a non-zero vector into a vector of unit length. */
    void normalize() {
        double l = magnitude();
        x /= l;
        y /= l;
        z /= l;
    }

    /* Multiplis this vector by the given scalar. */

    void operator*=(const double value) {

        x *= value;
        y *= value;
        z *= value;
    }

    /* Returns a copy of this vector scaled to the given value */

    Vector3 operator*(const double value) const {
        return Vector3(x * value, y * value, z * value);
    }

    /* Adds the given vector to this. */

    void operator+=(const Vector3 &v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    /* Returns the value of the given vector added to this. */
    Vector3 operator+(const Vector3 &v) const {

        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    bool operator<(const Vector3 &v) const {

        return ((x < v.x) && (y < v.y) && (z < v.z));
    }

    bool operator>(const Vector3 &v) const {

        return ((x > v.x) && (y > v.y) && (z > v.z));
    }

    /* Subtracts the given vector from this. */
    void operator-=(const Vector3 &v) {

        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    /* Returns the value of the given vector subtracted from this. */
    Vector3 operator-(const Vector3 &v) const {

        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    /* Adds the given vector to this, scaled by the given amount.
    We could use the preceding methods to accomplish this, but it is
    useful to have a separate method.*/
    void addScaledVector(const Vector3 &vector, double scale)

    {

        x += vector.x * scale;
        y += vector.y * scale;
        z += vector.z * scale;
    }

    /* Calculates and returns a component-wise product of this
            vector with the given vector.*/
    Vector3 componentProduct(const Vector3 &vector) const {

        return Vector3(x * vector.x, y * vector.y, z * vector.z);
    }

    /* Performs a component-wise product with the given vector
            and sets this vector to its result*/
    void componentProductUpdate(const Vector3 &vector) {

        x *= vector.x;
        y *= vector.y;
        z *= vector.z;
    }

    /* Calculates and returns the scalar product of this vector
            with the given vector*/

    double scalarProduct(const Vector3 &vector) const {

        return x * vector.x + y * vector.y + z * vector.z;
    }

    /* Calculates and returns the scalar product of this vector
    with the given vector. */

    double operator*(const Vector3 &vector) const {

        return x * vector.x + y * vector.y + z * vector.z;
    }

    /* Calculates and returns the vector product of this vector
            with the given vector. */

    Vector3 vectorProduct(const Vector3 &vector) const {

        return Vector3(y * vector.z - z * vector.y, z * vector.x - x * vector.z,
                       x * vector.y - y * vector.x);
    }

    /* Updates this vector to be the vector product of its
            current value and the given vector*/

    void operator%=(const Vector3 &vector) { *this = vectorProduct(vector); }

    /* Calculates and returns the vector product of this
            vector with the given vector. */

    Vector3 operator%(const Vector3 &vector) const {

        return Vector3(y * vector.z - z * vector.y, z * vector.x - x * vector.z,
                       x * vector.y - y * vector.x);
    }

    friend std::ostream &operator<<(std::ostream &out, Vector3 &v) {

        out << v.x << " " << v.y << " " << v.z;
        return out;
    }
};

class SquareCylinder_XY {

  private:
    // required info given to the constructor
    // region is not required
    double center_x;
    double center_y;
    double width;
    int region;
    double T2;
    double D;

    // calculated in constructor to save time in intersection algorithm
    double x_corner1;
    double y_corner1;
    double x_corner2;
    double y_corner2;

  public:
    SquareCylinder_XY(){};
    SquareCylinder_XY(double, double, double, int, double, double);
    ~SquareCylinder_XY(){};

    bool inside(Vector3 &);
    bool intersect(Line &, double &);
    Vector3 getNormal(double, double, double);
    int getregion(Vector3 &);
    double getT2(Vector3 &);
    double getD(Vector3 &);
};

class Sphere {

  private:
    double center_x;
    double center_y;
    double center_z;
    double radius;
    double T2;
    double D;
    int region;

  public:
    Sphere() {}
    Sphere(double, double, double, double, double, double, int);
    ~Sphere(){};

    bool inside(Vector3 &);
    bool intersect(Line &, double &);
    Vector3 getNormal(double, double, double);
    int getregion(Vector3 &);
    double getT2(Vector3 &);
    double getD(Vector3 &);
};

class Cylinder_XY {

  private:
    double center_x;
    double center_y;
    double radius;
    double D;
    double T2;
    int region;

  public:
    Cylinder_XY() {}
    Cylinder_XY(double, double, double, double, double, int);
    ~Cylinder_XY(){};

    bool inside(Vector3 &);
    bool intersect(Line &, double &);
    Vector3 getNormal(double, double, double);
    int getregion(Vector3 &);
    double getT2(Vector3 &);
    double getD(Vector3 &);
};

class MultiSphere {

  private:
    Sphere smallSphere;
    Sphere bigSphere;

  public:
    MultiSphere() {}
    MultiSphere(double, double, double, double, double, double, double, double,
                double, int, int);
    MultiSphere(Point &, double, double, double, double, double, double, int,
                int);
    ~MultiSphere(){};

    bool inside(Vector3 &);
    bool intersect(Line &, double &);
    Vector3 getNormal(double, double, double);
    int getregion(Vector3 &);
    double getT2(Vector3 &);
    double getD(Vector3 &);
};

class Line_XY {

  public:
    double _xi;
    double _yi;

    double _xf;
    double _yf;

    Line_XY(double, double, double, double);
    ~Line_XY(){};
    bool inside(Vector3 &);
    bool intersect(Line &, double &);
    bool intersect(Line_XY &, double &);
    Vector3 getNormal(double, double, double);
    int getregion(Vector3 &);
    double getT2(Vector3 &);
    double getD(Vector3 &);
};

#endif
