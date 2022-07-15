// We assume no holes -> a closed surface mesh

class Mesh {

  private:
    int size;
    Triangle list[MESH_NUM];

  public:
    __device__ __host__ Mesh() {}
    __device__ __host__ Mesh(int _size, real _T1, real _T2, real _D,
                             real _permeability, int _region) {
        size = _size;
        T1 = _T1;
        T2 = _T2;
        D = _D;
        permeability = _permeability;
        region = _region;
    }
    __device__ __host__ ~Mesh(){};

    __device__ __host__ bool inside(const Vector3 &r) const {

        Vector3 zhat(0.0, 0.0, 1.0);
        // projection on the xy plane
        Vector3 xyproj = r;
        xyproj.z = 0.0;
        Vector3 dr = zhat - xyproj;

        for (int i = 0; i < size; i++) {
            Vector3 pdr = list[i].v0 - xyproj;
            Vector3 intersection =
                xyproj + dr * ((normal * pdr) / (normal * dr));
            if (list[i].inside(intersection)) {
                noi++;
            }
        }
        return ((noi % 2) != 0);
    }

    __device__ __host__ bool intersect(const Vector3 &ri, const Vector3 &rf,
                                       real &v, Vector3 &n) const {

        double vBestEst = 10.0;
        double vtemp;

        for (int i = 0; i < size; i++) {
            if (list[i].intersect(ri, rf, vtemp) && vtemp < vBestEst) {
                vBestEst = vtemp;
                n = list[i].getNormal();
            }
        }

        v = vBestEst;
        return (vBestEst < 5.0);
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
