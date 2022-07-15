class Triangle {

  private:
    Vector3 p1;
    Vector3 p2;
    Vector3 p3;
    Vector3 normal;

  public:
    __device__ __host__ Triangle() {}
    __device__ __host__ Triangle(Vector3 _p1, Vector3 _p2, Vector3 _p3) {
        v0 = _p1;
        v1 = _p2;
        v2 = _p3;
        normal = (v1 - v0) % (v2 - v0);
        normal.normalize();
    }

    __device__ __host__ ~Triangle() {}

    __device__ __host__ bool intersect(const Vector3 &ri, const Vector3 &rf,
                                       real &v) const {
        // intersection with entire plane defined by triangle
        Vector3 dr = rf - ri;
        Vector3 pdr = v0 - ri;
        real v = (normal * pdr) / (normal * dr); // intersection point
        if (v < 0.0 || v > 1.0 || doub_equal(v * dr.magnitude(), 0.0)) {
            return false;
        };

        // intersection point
        Vector3 intP = ri + dr * v;
        return inside(intP);
    }

    __device__ __host__ bool inside(const Vector3 &r) {

        Vector3 w = r - v0;
        Vector3 u = v1 - v0;
        Vector3 v = v2 - v0;

        real udotv = u * v;
        real wdotv = w * v;
        real wdotu = w * u;
        real vmag = v * v;
        real umag = u * u;
        real denom = udotv * udotv - umag * vmag;
        real s1 = (udotv * wdotv - vmag * wdotu) / denom;
        real t1 = (udotv * wdotu - umag * wdotv) / denom;

        if (s1 > 0.0 && t1 > 0.0 && (s1 + t1) < 1.0) {
            return true;
        }
        return false;
    }

    // here r is a point on the surface
    __device__ __host__ Vector3 getNormal(Vector3 &r) const { return normal; }
};
