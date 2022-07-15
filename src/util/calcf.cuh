template <class Basis>
__global__ void kernelCalcf(curandState *globalState, Basis *basis,
                            Lattice *lat, int *insideBasis, int points) {

    const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
    curandState localState = globalState[tid];
    Vector3 r;
    lat->initializeUniformly(r, localState);
    insideBasis[tid] = (lat->inRegion(basis, r) > 0);
}

template <class Basis>
__host__ double calcf(Basis *basis, Lattice &lattice, int points, int seed,
                      int blocks, int threads, int numOfSM = 14,
                      int maxThreads = 1024) {

    cudaVector<curandState> devStates(points);
    cudaVector<Basis> dev_basis(lattice.getBasisSize());
    cudaScalar<Lattice> dev_lat(lattice);
    cudaVector<int> insideBasis(points);
    cudaScalar<int> sum;
    sum.malloc();

    for (int i = 0; i < lattice.getBasisSize(); i++) {
        dev_basis[i] = basis[i];
    }

    setup_kernel<<<blocks, threads>>>(devStates.getPointer(), seed);

    dev_basis.copyToDevice();
    dev_lat.copyToDevice();

    kernelCalcf<Basis><<<blocks, threads>>>(
        devStates.getPointer(), dev_basis.getPointer(), dev_lat.getPointer(),
        insideBasis.getPointer(), points);

    insideBasis.sum(sum, maxthreads, numOfSM, 0);
    return sum.getValue() / points;
}