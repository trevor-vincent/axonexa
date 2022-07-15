// mainly for testing purposes
class Empty {

  private:
  public:
    __device__ __host__ Empty() {}
    __device__ __host__ ~Empty(){};

    __device__ __host__ bool inside(const Vector3 &r) const { return false; }

    __device__ __host__ bool intersect(const Vector3 ri, const Vector3 rf,
                                       real *v) const {
        return false;
    }

    // here r is a point on the surface
    __device__ __host__ Vector3 getNormal(Vector3 &r) const {
        return Vector3(0.0, 0.0, 0.0);
    }

    __device__ __host__ int getRegion(const Vector3 &r) const { return 0; }

    __device__ __host__ real getT2(const Vector3 &r) const { return -1.0; }

    __device__ __host__ real getD(const Vector3 &r) const { return -1.0; }

    __device__ __host__ real getPermeability() const { return 0.0; }
};
