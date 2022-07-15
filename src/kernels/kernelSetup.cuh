template <class Params>
__global__ void setup_kernel(curandState *state, const Params *par) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(par->seed, tid, 0, &state[tid]);
}

__global__ void setup_kernel(curandState *state, int seed) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, tid, 0, &state[tid]);
}

__device__ int randomSign(curandState &localState) {
    return (curand(&localState) % 2) ? 1 : -1;
}

__host__ int randomSign() { return (rand() % 2) ? 1 : -1; }