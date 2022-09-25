#pragma once

__device__ double atomicAddInHouse(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(
            address_as_ull, assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

// number of blocks should equal the number of SMs
// number of threads should equal the maximum allowable number of threads in a
// block
template <class T>
__global__ void _functionReduceAtom(T *sum, T *g_idata, const unsigned int n) {

    T myVal = *sum = (T)0;

    const int numThreads =
        blockDim.x; // should = maximum number of threads allowed on SM

    { // 1) Use fastest memory first.
        const int gridSize = numThreads * gridDim.x;

        // for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i +=
        // gridSize) myVal = fcn1(fcn(i), myVal);
        for (int i = n - 1 - (blockIdx.x * blockDim.x + threadIdx.x); i >= 0;
             i -= gridSize)
            myVal += g_idata[i];
    }

    // 2) Use the second fastest memory (shared memory) in a warp
    // synchronous fashion.
    // Create shared memory for per-block reduction.
    // Reuse the registers in the first warp.
    volatile extern __shared__ T smem[];

    // put all the register values into a shared memory
    if (threadIdx.x >= WARP_SIZE)
        smem[threadIdx.x - WARP_SIZE] = myVal;
    __syncthreads(); // wait for all threads in the block to complete.

    if (threadIdx.x < WARP_SIZE) {
        // now using just one warp. The SM can only run one warp at a time

#pragma unroll
        for (int i = threadIdx.x; i < (numThreads - WARP_SIZE); i += WARP_SIZE)
            myVal += (T)smem[i];
        smem[threadIdx.x] =
            myVal; // save myVal in this warp to the start of smem
    }

    // reduce shared memory.
    if (threadIdx.x < 16)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 16];
    if (threadIdx.x < 8)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 8];
    if (threadIdx.x < 4)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 4];
    if (threadIdx.x < 2)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 2];
    if (threadIdx.x < 1)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 1];

    // 3) Use global memory as a last resort to transfer results to the host
    // write result for each block to global mem

    if (threadIdx.x == 0)
        atomicAddInHouse(sum, smem[0]);
}

template <class T, class Transform>
__global__ void _functionTransformAndReduceAtom(T *sum, T *g_idata,
                                                const unsigned int n) {

    T myVal = *sum = (T)0;
    Transform trans;

    const int numThreads =
        blockDim.x; // should = maximum number of threads allowed on SM

    { // 1) Use fastest memory first.
        const int gridSize = numThreads * gridDim.x;

        // for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i +=
        // gridSize) myVal = fcn1(fcn(i), myVal);
        for (int i = n - 1 - (blockIdx.x * blockDim.x + threadIdx.x); i >= 0;
             i -= gridSize)
            myVal += trans(g_idata[i]);
    }

    // 2) Use the second fastest memory (shared memory) in a warp
    // synchronous fashion.
    // Create shared memory for per-block reduction.
    // Reuse the registers in the first warp.
    volatile extern __shared__ T smem[];

    // put all the register values into a shared memory
    if (threadIdx.x >= WARP_SIZE)
        smem[threadIdx.x - WARP_SIZE] = myVal;
    __syncthreads(); // wait for all threads in the block to complete.

    if (threadIdx.x < WARP_SIZE) {
        // now using just one warp. The SM can only run one warp at a time

#pragma unroll
        for (int i = threadIdx.x; i < (numThreads - WARP_SIZE); i += WARP_SIZE)
            myVal += (T)smem[i];
        smem[threadIdx.x] =
            myVal; // save myVal in this warp to the start of smem
    }

    // reduce shared memory.
    if (threadIdx.x < 16)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 16];
    if (threadIdx.x < 8)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 8];
    if (threadIdx.x < 4)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 4];
    if (threadIdx.x < 2)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 2];
    if (threadIdx.x < 1)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 1];

    // 3) Use global memory as a last resort to transfer results to the host
    // write result for each block to global mem

    if (threadIdx.x == 0)
        atomicAddInHouse(sum, smem[0]);
}

template <class T, class Transform>
__global__ void _functionTransformAndSumTwoVectorsAtom(T *sum, T *g_idata1,
                                                       T *g_idata2,
                                                       const unsigned int n) {

    T myVal = *sum = (T)0;
    Transform trans;

    const int numThreads =
        blockDim.x; // should = maximum number of threads allowed on SM

    { // 1) Use fastest memory first.
        const int gridSize = numThreads * gridDim.x;

        // for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i +=
        // gridSize) myVal = fcn1(fcn(i), myVal);
        for (int i = n - 1 - (blockIdx.x * blockDim.x + threadIdx.x); i >= 0;
             i -= gridSize)
            myVal += trans(g_idata1[i], g_idata2[i]);
    }

    // 2) Use the second fastest memory (shared memory) in a warp
    // synchronous fashion.
    // Create shared memory for per-block reduction.
    // Reuse the registers in the first warp.
    volatile extern __shared__ T smem[];

    // put all the register values into a shared memory
    if (threadIdx.x >= WARP_SIZE)
        smem[threadIdx.x - WARP_SIZE] = myVal;
    __syncthreads(); // wait for all threads in the block to complete.

    if (threadIdx.x < WARP_SIZE) {
        // now using just one warp. The SM can only run one warp at a time

#pragma unroll
        for (int i = threadIdx.x; i < (numThreads - WARP_SIZE); i += WARP_SIZE)
            myVal += (T)smem[i];
        smem[threadIdx.x] =
            myVal; // save myVal in this warp to the start of smem
    }

    // reduce shared memory.
    if (threadIdx.x < 16)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 16];
    if (threadIdx.x < 8)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 8];
    if (threadIdx.x < 4)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 4];
    if (threadIdx.x < 2)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 2];
    if (threadIdx.x < 1)
        smem[threadIdx.x] += (T)smem[threadIdx.x + 1];

    // 3) Use global memory as a last resort to transfer results to the host
    // write result for each block to global mem

    if (threadIdx.x == 0)
        atomicAddInHouse(sum, smem[0]);
}
