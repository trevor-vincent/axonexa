class Sphere {

  private:
    real centerX;
    real centerY;
    real centerZ;
    real radius;
    real D;
    real T2;
    real T1;
    real permeability;
    real sphereEPS;
    int region;

  public:
    __device__ __host__ Sphere() {}
    __device__ __host__ Sphere(real _centerX, real _centerY, real _centerZ,
                               real _radius, real _T2, real _T1, real _D,
                               int _region, real _permeability) {

        centerX = _centerX;
        centerY = _centerY;
        centerZ = _centerZ;
        radius = _radius;
        D = _D;
        T2 = _T2;
        T1 = _T1;
        region = _region;
        permeability = _permeability;
        sphereEPS = 1E-14;
    }

    __device__ __host__ Sphere(real _centerX, real _centerY, real _centerZ,
                               real _radius, real _T2, real _T1, real _D,
                               int _region, real _permeability, real _eps) {

        centerX = _centerX;
        centerY = _centerY;
        centerZ = _centerZ;
        radius = _radius;
        D = _D;
        T2 = _T2;
        T1 = _T1;
        region = _region;
        permeability = _permeability;
        sphereEPS = _eps;
    }

    __device__ __host__ ~Sphere(){};

    // this will need to be changed once we go to full generality
    __device__ Vector3 unifRand(curandState localState) const {

        real x = radius * (2.0 * curand_uniform(&localState) - 1);
        real y = radius * (2.0 * curand_uniform(&localState) - 1);
        real z = radius * (2.0 * curand_uniform(&localState) - 1);
        Vector3 r(x, y, z);

        while (!inside(r)) {

            r = Vector3(radius * (2.0 * curand_uniform(&localState) - 1),
                        radius * (2.0 * curand_uniform(&localState) - 1),
                        radius * (2.0 * curand_uniform(&localState) - 1));
        }

        return r;
    }

    __host__ Vector3 unifRandCPU() const {

        real x = radius * (2.0 * unifRandCPP() - 1);
        real y = radius * (2.0 * unifRandCPP() - 1);
        real z = radius * (2.0 * unifRandCPP() - 1);
        Vector3 r(x, y, z);

        while (!inside(r)) {

            r = Vector3(radius * (2.0 * unifRandCPP() - 1),
                        radius * (2.0 * unifRandCPP() - 1),
                        radius * (2.0 * unifRandCPP() - 1));
        }

        return r;
    }

    __device__ __host__ bool inside(const Vector3 &r) const {
        return ((r.x - centerX) * (r.x - centerX) +
                    (r.y - centerY) * (r.y - centerY) +
                    (r.z - centerZ) * (r.z - centerZ) <
                radius * radius);
    }

    __device__ __host__ bool inside(real x, real y, real z) const {
        return ((x - centerX) * (x - centerX) + (y - centerY) * (y - centerY) +
                    (z - centerZ) * (z - centerZ) <
                radius * radius);
    }

    __device__ __host__ bool intersect(const Vector3 &ri, const Vector3 &rf,
                                       real &v) const {

        Vector3 dr = rf - ri;
        real step_mag = dr.magnitude();

        real a = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
        real b = 2.0 * ri.x * dr.x - 2.0 * dr.x * centerX + 2.0 * ri.y * dr.y -
                 2.0 * dr.y * centerY + 2.0 * ri.z * dr.z -
                 2.0 * dr.z * centerZ;
        real c = ri.x * ri.x + ri.y * ri.y + ri.z * ri.z - 2 * ri.x * centerX -
                 2 * ri.y * centerY - 2 * ri.z * centerZ + centerX * centerX +
                 centerY * centerY + centerZ * centerZ - radius * radius;

        real q = -.5 * (b + sgn(b) * sqrt(b * b - 4 * a * c));
        real root1 = q / a;
        real root2 = c / q;

        bool s1 = (root1 > 0.0 && root1 < 1.0 && b * b > 4 * a * c &&
                   !real_equal(root1 * step_mag, 0.0, sphereEPS));
        bool s2 = (root2 > 0.0 && root2 < 1.0 && b * b > 4 * a * c &&
                   !real_equal(root2 * step_mag, 0.0, sphereEPS));
        bool s3 = (fabs(root1) < fabs(root2));

        if ((s1 && s2 && s3) || (s1 && !s2)) {
            v = root1;
            return true;
        }

        else if ((s1 && s2 && !s3) || (s2 && !s1)) {
            v = root2;
            return true;
        }

        else {
            return false;
        }
    }

    // here r is a point on the surface
    __device__ __host__ Vector3 getNormal(Vector3 &r) const {

        double n_x = r.x - centerX;
        double n_y = r.y - centerY;
        double n_z = r.z - centerZ;
        double mag = sqrt(n_x * n_x + n_y * n_y + n_z * n_z);
        return Vector3(n_x / mag, n_y / mag, n_z / mag);
    }

    __device__ __host__ real getRadius() const { return radius; }

    __device__ __host__ int getRegion(const Vector3 &r) const { return region; }

    __device__ __host__ real getT2(const Vector3 &r) const {
        if (inside(r)) {
            return T2;
        }
        return -1.0;
    }

    __device__ __host__ real getD(const Vector3 &r) const {
        if (inside(r)) {
            return D;
        }
        return -1.0;
    }

    __device__ __host__ real getD() const { return D; }

    __device__ __host__ real getPermeability() const { return permeability; }

    __device__ __host__ Vector3 getCenter() const {
        return Vector3(centerX, centerY, centerZ);
    }

    __host__ void setCenter(Vector3 v) {
        centerX = v.x;
        centerY = v.y;
        centerZ = v.z;
    }

    __host__ void setRadius(real _r) { radius = _r; }

    __host__ void setEPS(real _sphereEPS) { sphereEPS = _sphereEPS; }

    __host__ void setRegion(int _region) { region = _region; }

    __host__ int getRegion() { return region; }

    __device__ void randUnif(Vector3 &r, curandState &localState) const {
        do {
            r = Vector3((2.0 * curand_uniform(&localState) - 1.0) * radius,
                        (2.0 * curand_uniform(&localState) - 1.0) * radius,
                        (2.0 * curand_uniform(&localState) - 1.0) * radius) +
                getCenter();
        } while (!inside(r));
    }

    __host__ void randUnif(Vector3 &r) const {
        do {
            r = Vector3((2.0 * unifRandCPP() - 1.0) * radius,
                        (2.0 * unifRandCPP() - 1.0) * radius,
                        (2.0 * unifRandCPP() - 1.0) * radius) +
                getCenter();
        } while (!inside(r));
    }
};
