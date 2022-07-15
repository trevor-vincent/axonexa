// We assume no holes -> a closed surface DamagedAxon

class DamagedAxon {

  private:
    Cylinder_XY cyl;
    Sphere spher;

  public:
    __device__ __host__ DamagedAxon() {}
    __device__ __host__ DamagedAxon(Vector3 _center, real _radiusCyl,
                                    real _radiusSphere, real _T1, real _T2,
                                    real _D, real _permeability, int _region) {
        cyl = Cylinder_XY(_center.x, _center.y, _radiusCyl, _T2, _T1, _D,
                          _region, _permeability);
        spher = Sphere(_center.x, _center.y, _center.z, _radiusSphere, _T2, _T1,
                       _D, _region, _permeability);
    }
    __device__ __host__ ~DamagedAxon(){};

    __device__ __host__ bool inside(const Vector3 &r) const {
        return (cyl.inside(r) || spher.inside(r));
    }

    __device__ __host__ bool intersect(const Vector3 &ri, const Vector3 &rf,
                                       real &v, Vector3 &n) const {

        real v1, v2;
        Vector3 n1, n2;

        bool s1 = cyl.intersect(ri, rf, v1, n1);
        bool s2 spher.intersect(ri, rf, v2, n2);
        if (!s1 && !s2) {
            return false;
        }

        Vector3 rnew1 = (rf - ri) * v1 + ri;
        Vector3 rnew2 = (rf - ri) * v2 + ri;

        bool s3 = v1 < v2;
        bool s4 = spher.inside(rnew1);
        bool s5 = cyl.inside(rnew2);
        bool s6 = inside(ri);

        // cylinder intersection but inside spher
        if (s1 && (!s2 || s3)) {
            if (s4 && s6) {
                return false;
            } else {
                v = v1;
                n = n1;
                return true;
            }
        }

        else if (s2 && (!s1 || !s3)) {
            if (s5 && s6) {
                return false;
            } else {
                v = v2;
                n = n2;
                return true;
            }
        }
    }

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
};
