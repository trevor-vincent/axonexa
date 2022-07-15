class Vector3 {

  public:
    /* Holds the value along the x axis. */
    real x;

    /* Holds the value along the y axis */
    real y;

    /* Holds the value along the z axis */
    real z;

    /* The default constructor creates a zero vector */
    __device__ __host__ Vector3() : x(0), y(0), z(0) {}

    /* The explicit constructor creates a vector with the
    given components		*/

    __device__ __host__ Vector3(const real x, const real y, const real z)
        : x(x), y(y), z(z) {}

    // /* Flips all the components of the vector */

    __device__ __host__ void invert() {

        x = -x;
        y = -y;
        z = -z;
    }

    /* Gets the magnitude of this vector. */
    __device__ __host__ real magnitude() const {

        return sqrt(x * x + y * y + z * z);
    }

    /* Gets the squared magnitude of this vector. In some
    cases we do not need the exact magnitude; for example, when
    we need to compare two magnitudes to see which is greater, it
    is faster the compare the squares of the magnitudes. This is
    b/c we do not need to call the sqrt function.*/
    __device__ __host__ real squareMagnitude() const {

        return x * x + y * y + z * z;
    }

    /* Turns a non-zero vector into a vector of unit length. */
    __device__ __host__ void normalize() {
        real l = magnitude();
        x /= l;
        y /= l;
        z /= l;
    }

    /* Multiplis this vector by the given scalar. */

    __device__ __host__ void operator*=(const real value) {

        x *= value;
        y *= value;
        z *= value;
    }

    /* Returns a copy of this vector scaled to the given value */

    __device__ __host__ Vector3 operator*(const real value) const {
        return Vector3(x * value, y * value, z * value);
    }

    __device__ __host__ Vector3 operator/(const real value) const {
        return Vector3(x / value, y / value, z / value);
    }

    /* Adds the given vector to this. */

    __device__ __host__ void operator+=(const Vector3 &v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    /* Returns the value of the given vector added to this. */
    __device__ __host__ Vector3 operator+(const Vector3 &v) const {

        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    __device__ __host__ bool operator<(const Vector3 &v) const {

        return ((x < v.x) && (y < v.y) && (z < v.z));
    }

    __device__ __host__ bool operator>(const Vector3 &v) const {

        return ((x > v.x) && (y > v.y) && (z > v.z));
    }

    /* Subtracts the given vector from this. */
    __device__ __host__ void operator-=(const Vector3 &v) {

        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    /* Returns the value of the given vector subtracted from this. */
    __device__ __host__ Vector3 operator-(const Vector3 &v) const {

        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    /* Adds the given vector to this, scaled by the given amount.
    We could use the preceding methods to accomplish this, but it is
    useful to have a separate method.*/
    __device__ __host__ void addScaledVector(const Vector3 &vector, real scale)

    {

        x += vector.x * scale;
        y += vector.y * scale;
        z += vector.z * scale;
    }

    /* Calculates and returns a component-wise product of this
            vector with the given vector.*/
    __device__ __host__ Vector3 componentProduct(const Vector3 &vector) const {

        return Vector3(x * vector.x, y * vector.y, z * vector.z);
    }

    /* Performs a component-wise product with the given vector
            and sets this vector to its result*/
    __device__ __host__ void componentProductUpdate(const Vector3 &vector) {

        x *= vector.x;
        y *= vector.y;
        z *= vector.z;
    }

    /* Calculates and returns the scalar product of this vector
            with the given vector*/

    __device__ __host__ real scalarProduct(const Vector3 &vector) const {

        return x * vector.x + y * vector.y + z * vector.z;
    }

    /* Calculates and returns the scalar product of this vector
    with the given vector. */

    __device__ __host__ real operator*(const Vector3 &vector) const {

        return x * vector.x + y * vector.y + z * vector.z;
    }

    /* Calculates and returns the vector product of this vector
            with the given vector. */

    __device__ __host__ Vector3 vectorProduct(const Vector3 &vector) const {

        return Vector3(y * vector.z - z * vector.y, z * vector.x - x * vector.z,
                       x * vector.y - y * vector.x);
    }

    /* Updates this vector to be the vector product of its
            current value and the given vector*/

    __device__ __host__ void operator%=(const Vector3 &vector) {

        *this = vectorProduct(vector);
    }

    /* Calculates and returns the vector product of this
            vector with the given vector. */

    __device__ __host__ Vector3 operator%(const Vector3 &vector) const {

        return Vector3(y * vector.z - z * vector.y, z * vector.x - x * vector.z,
                       x * vector.y - y * vector.x);
    }

    __host__ friend std::ostream &operator<<(std::ostream &out, Vector3 &v) {

        out << v.x << " " << v.y << " " << v.z;
        return out;
    }
};
