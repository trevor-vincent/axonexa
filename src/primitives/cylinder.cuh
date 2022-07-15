class Cylinder {

  private:
    Vector3 center_cyl; // In Simulation Coordinates
    Vector3 ahat;       // cylinder x axis in Simulation Coordinates
    Vector3 bhat;       // cylinder y axis in Simulation Coordinates
    Vector3 chat;       // cylinder axis in Simulation Coordinates
    real radius;
    real T2;
    real T1;
    real D;
    real permeability

        public :

        __device__ __host__
        Cylinder() {}
    __device__ __host__ Cylinder(Vector3 _center_simul, Vector3 _axis,
                                 real _radius, real _T2, real _T1, real _D,
                                 int _region, real _permeability) {

        chat = _axis;
        chat.normalize();

        // cross the axis into the z-axis to get the cylinder x-axis
        ahat = chat % Vector3(0.0, 0.0, 1.0);

        if (doub_equal(ahat.x, 0.0) && doub_equal(ahat.y, 0.0) &&
            doub_equal(ahat.z, 0.0)) {
            // you should use Cylinder_XY if this happens
            ahat = Vector3(1.0, 0.0, 0.0);
            bhat = Vector3(0.0, 1.0, 0.0);
        }

        else {
            ahat.normalize();
            // cross the cylinder axis into the x-axis to get the y-axis
            bhat = chat % ahat;
            bhat.normalize;
        }

        radius = _radius;
        T2 = _T2;
        T1 = _T1;
        D = _D;
        region = _region;
        permeability = _permeability;
        center_cyl = getCylinderCoordinates(
            _center_simul); // center of cylinder in cylinder coordinates
                            // (ahat,bhat,chat basis)
    }

    __device__ __host__ ~Cylinder(){};

    __device__ __host__ Vector3 getCylinderCoordinates(Vector3 &r_simul) {

        // r_simul in cylinder coordinates, i.e, in the ahat,bhat,chat basis
        Vector3 r_cyl;
        r_cyl.x = r_simul * ahat;
        r_cyl.y = r_simul * bhat;
        r_cul.z = r_simul * chat;
        return r_cyl;
    }

    __device__ __host__ Vector3 getSimulationCoordinates(Vector3 &r_cyl) {

        // r_cyl in simulation coordinates, i.e, in the xhat,yhat,zhat basis
        Vector3 r_simul;
        r_simul.x = r_cyl * Vector3(1.0, 0.0, 0.0);
        r_simul.y = r_cyl * Vector3(0.0, 1.0, 0.0);
        r_simul.z = r_cyl * Vector3(0.0, 0.0, 1.0);
        return r_simul;
    }

    __device__ Vector3 unifRand(curandState localState) const {

        Vector3 r_simul, r_cyl;

        do {
            r_cyl.x =
                radius * (2.0 * curand_uniform(&localState) - 1) + center_cyl.x;
            r_cyl.y =
                radius * (2.0 * curand_uniform(&localState) - 1) + center_cyl.y;
            r_cyl.z = 2.0 * radius * curand_uniform(&localState) + center_cyl.z;
            r_simul = getSimulCoordinates(r_cyl);
        } while (!inside(r_simul));

        return r;
    }

    __host__ Vector3 unifRandCPU() const {

        Vector3 r_simul, r_cyl;

        do {
            r_cyl.x = radius * (2.0 * unifRandCPP() - 1) + center_cyl.x;
            r_cyl.y = radius * (2.0 * unifRandCPP() - 1) + center_cyl.y;
            r_cyl.z = 2.0 * radius * unifRandCPP() + center_cyl.z;
            r_simul = getSimulCoordinates(r_cyl);
        } while (!inside(r_simul));

        return r_simul;
    }

    __device__ __host__ bool inside(Vector3 &r_simul) const {

        Vector3 r_cyl = getCylinderCoordinates(r_simul);
        real a = r_cyl.x - center_cyl.x;
        real b = r_cyl.y - center_cyl.y;
        return (a * a + b * b < radius * radius);
    }

    __device__ __host__ bool intersect(const Vector3 ri_simul,
                                       const Vector3 rf_simul, real *v) const {

        Vector3 ri = getCylinderCoordinates(ri_simul);
        Vector3 rf = getCylinderCoordinates(rf_simul);

        Vector3 dr = rf - ri;
        real step_mag = dr.magnitude();

        // real a = dr.x*dr.x + dr.y*dr.y;
        // real b = 2*ri.x*dr.x + 2*ri.y*dr.y;
        // real c = ri.x*ri.x + ri.y*ri.y - radius*radius;

        real a = dr.x * dr.x + dr.y * dr.y;
        real b = 2.0 * ri.x * dr.x - 2.0 * dr.x * center_cyl.x +
                 2.0 * ri.y * dr.y - 2.0 * dr.y * center_cyl.y;
        real c = ri.x * ri.x + ri.y * ri.y - 2.0 * ri.x * center_cyl.x -
                 2.0 * ri.y * center_cyl.y + center_cyl.x * center_cyl.x +
                 center_cyl.y * center_cyl.y - radius * radius;

        real q = -.5 * (b + sgn(b) * sqrt(b * b - 4 * a * c));
        real root1 = q / a;
        real root2 = c / q;

        bool s1 = (root1 > 0.0 && root1 < 1.0 && b * b > 4 * a * c &&
                   !doub_equal(root1 * step_mag, 0.0));
        bool s2 = (root2 > 0.0 && root2 < 1.0 && b * b > 4 * a * c &&
                   !doub_equal(root2 * step_mag, 0.0));
        bool s3 = (fabs(root1) < fabs(root2));

        if ((s1 && s2 && s3) || (s1 && !s2)) {
            *v = root1;
            return true;
        }

        else if ((s1 && s2 && !s3) || (s2 && !s1)) {
            *v = root2;
            return true;
        }

        else {
            return false;
        }
    }

    // here r is a point on the surface
    __device__ __host__ Vector3 getNormal(Vector3 &r_simul) const {
        Vector3 r_cyl = getCylinderCoordinates(r_simul);
        Vector3 n_cyl, n_simul;
        n_cyl.x = r_cyl.x - center_cyl.x;
        n_cyl.y = r_cyl.y - center_cyl.y;
        n_cyl.z = 0.0;
        n_simul = getSimulCoordinates(n_cyl);
        n_simul.normalize();
        return n_simul;
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

    __device__ __host__ real getPermeability() const { return permeability; }

    __device__ __host__ Vector3 getCenter() {
        return getSimulCoordinates(center_cyl);
    }

    __host__ void setCenter(Vector3 &v) {
        center_cyl = getCylinderCoordinates(v);
    }

    __host__ void setRadius(real _r) { radius = _r; }
};
