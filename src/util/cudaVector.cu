#include "reduction.cuh"
#include "cudaVector.cuh"

// log all errors into logger

template <class T> cudaVector<T>::cudaVector() {
    alloc = false;
    copied = false;
    device = 0;
}

template <class T> cudaVector<T>::cudaVector(int size) {
    _v.resize(size);
    safe_cuda(cudaMalloc((void **)&dev_ptr, size * sizeof(T)));
    alloc = true;
    copied = false;
    device = 0;
}

template <class T> cudaVector<T>::cudaVector(int size, T val) {
    _v = std::vector<T>(size, val);
    safe_cuda(cudaMalloc((void **)&dev_ptr, size * sizeof(T)));
    alloc = true;
    copied = false;
    device = 0;
}

template <class T> cudaVector<T>::~cudaVector() {

    if (alloc) {
        safe_cuda(cudaSetDevice(device));
        safe_cuda(cudaFree(dev_ptr));
    }
}

template <class T> void cudaVector<T>::setDevice(int _device) {
    device = _device;
}

template <class T> void cudaVector<T>::copyToDevice() {
    if (alloc == false) {
        safe_cuda(cudaMalloc((void **)&dev_ptr, _v.size() * sizeof(T)));
        alloc = true;
    }
    safe_cuda(cudaMemcpy(dev_ptr, &_v[0], _v.size() * sizeof(T),
                         cudaMemcpyHostToDevice));
    copied = true;
}

template <class T> void cudaVector<T>::copyFromDevice() {
    if (alloc == true) {
        safe_cuda(cudaMemcpy(&_v[0], dev_ptr, _v.size() * sizeof(T),
                             cudaMemcpyDeviceToHost));
    }
}

template <class T> void cudaVector<T>::malloc(int asize) {
    safe_cuda(cudaMalloc((void **)&dev_ptr, asize * sizeof(T)));
    alloc = true;
    _v.resize(asize);
}

template <class T> int cudaVector<T>::size() { return _v.size(); }

template <class T> void cudaVector<T>::resize(int size) { _v.resize(size); }

template <class T> T &cudaVector<T>::operator[](int i) { return _v[i]; }

template <class T> void cudaVector<T>::operator=(pinnedVector<T> &v) {

    int iter;
    if (v.size() >= _v.size()) {
        iter = _v.size();
    } else {
        iter = v.size();
    }

    for (int i = 0; i < iter; i++) {
        _v[i] = v[i];
    }
}
template <class T> T *cudaVector<T>::getPointer() { return dev_ptr; }

template <class T> void cudaVector<T>::copyTo(std::vector<T> &v) {

    int iter;
    if (v.size() >= _v.size()) {
        iter = _v.size();
    } else {
        iter = v.size();
    }

    for (int i = 0; i < iter; i++) {
        v[i] = _v[i];
    }
}

template <class T> void cudaVector<T>::operator=(std::vector<T> &v) {

    int iter;
    if (v.size() >= _v.size()) {
        iter = _v.size();
    } else {
        iter = v.size();
    }

    for (int i = 0; i < iter; i++) {
        _v[i] = v[i];
    }
}

template <class T>
void cudaVector<T>::sum(cudaScalar<T> &sum, int threads, int blocks,
                        cudaStream_t stream) {
    _functionReduceAtom<<<blocks, threads, sizeof(T) * (threads - WARP_SIZE),
                          stream>>>(
        sum.getPointer(), cudaVector<T>::getPointer(), cudaVector<T>::size());
}

// sum ( X )
template <class T>
void cudaVector<T>::sum(cudaVector<T> &sum, int threads, int blocks, int chunks,
                        int sizeOfChunk, cudaStream_t stream) {
    for (int i = 0; i < chunks; i++) {
        _functionReduceAtom<<<blocks, threads,
                              sizeof(T) * (threads - WARP_SIZE), stream>>>(
            sum.getPointer() + i, cudaVector<T>::getPointer() + i * sizeOfChunk,
            sizeOfChunk);
    }
}

// sum ( transform ( X ) )
template <class T>
template <class Transform>
void cudaVector<T>::transformAndSum(cudaScalar<T> &sum, int threads, int blocks,
                                    cudaStream_t stream) {
    _functionTransformAndReduceAtom<T, Transform>
        <<<blocks, threads, sizeof(T) * (threads - WARP_SIZE), stream>>>(
            sum.getPointer(), cudaVector<T>::getPointer(),
            cudaVector<T>::size());
}

// sum ( transform ( X ) )
template <class T>
template <class Transform>
void cudaVector<T>::transformAndSum(cudaVector<T> &sum, int threads, int blocks,
                                    int chunks, int sizeOfChunk,
                                    cudaStream_t stream) {
    for (int i = 0; i < chunks; i++) {
        _functionTransformAndReduceAtom<T, Transform>
            <<<blocks, threads, sizeof(T) * (threads - WARP_SIZE), stream>>>(
                sum.getPointer() + i,
                cudaVector<T>::getPointer() + i * sizeOfChunk, sizeOfChunk);
    }
}

// sum ( transform ( X , Y ) )
template <class T>
template <class Transform>
void cudaVector<T>::transformAndSumTwoVectors(cudaVector<T> &sum,
                                              cudaVector<T> &secondVector,
                                              int threads, int blocks,
                                              int chunks, int sizeOfChunk,
                                              cudaStream_t stream) {
    for (int i = 0; i < chunks; i++) {
        _functionTransformAndSumTwoVectorsAtom<T, Transform>
            <<<blocks, threads, sizeof(T) * (threads - WARP_SIZE), stream>>>(
                sum.getPointer() + i,
                cudaVector<T>::getPointer() + i * sizeOfChunk,
                secondVector.getPointer(), sizeOfChunk);
    }
}
