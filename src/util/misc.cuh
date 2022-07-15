__device__ __host__ int sgn(real val) { return (0.0 < val) - (val < 0.0); }

inline void safe_cuda(cudaError_t a) {
    if (a != cudaSuccess) {
        printf(" \n CUDA Error: %s (err_num=%d) \n", cudaGetErrorString(a), a);
        cudaDeviceReset();
        assert(0);
    }
}

void printDeviceInfo(int device_num) {

    struct cudaDeviceProp device_prop;
    safe_cuda(cudaGetDeviceProperties(&device_prop, device_num));
    std::cout << " You are using Device " << device_num << " : "
              << device_prop.name << std::endl;
}

template <typename T> struct COS_UNARY {
    __host__ __device__ T operator()(const T &x) const { return cos(x); }
};

template <typename T> struct SIN_UNARY {
    __host__ __device__ T operator()(const T &x) const { return sin(x); }
};

template <typename T> struct COS_EXP_BINARY {
    __host__ __device__ T operator()(const T &x1, const T &x2) const {
        return cos(x1) * exp(x2);
    }
};

template <typename T> struct SIN_EXP_BINARY {
    __host__ __device__ T operator()(const T &x1, const T &x2) const {
        return sin(x1) * exp(x2);
    }
};

double linear_regression(std::vector<real> y, std::vector<real> x) {

    if (y.size() != x.size()) {
        std::cout << "ERROR LINEAR_REGRESSION Y.SIZE != X.SIZE " << std::endl;
        exit(0);
    }
    double term1 = y.size(), term2 = 0.0, term3 = 0.0, term4 = 0.0, term5 = 0.0,
           delta;

    for (int i = 0; i < y.size(); i++) {
        term2 = term2 + (x[i] * y[i]);
        term3 = term3 + (x[i]);
        term4 = term4 + (y[i]);
        term5 = term5 + (x[i] * x[i]);
    }

    delta = term1 * term5 - term3 * term3;
    return (1 / delta) * (term1 * term2 - term3 * term4);
}

// only call this function to debug. It forces synchronization of the CPU and
// the GPU
__host__ void cudaLastErrorCheck() {

    cudaDeviceSynchronize();
    if (cudaPeekAtLastError() != cudaSuccess) {

        std::cout << std::endl
                  << " Last Error: " << std::endl
                  << std::endl
                  << cudaGetErrorString(cudaGetLastError()) << std::endl;
        cudaDeviceReset();
        exit(0);
    }
}

//
// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRandCPP() { return rand() / double(RAND_MAX); }

//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRandCPP(double a, double b) { return (b - a) * unifRandCPP() + a; }
