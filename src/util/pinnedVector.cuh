#ifndef PINNEDVECTOR_H_
#define PINNEDVECTOR_H_

template <class X> class cudaVector;

template <class T> class pinnedVector {

  private:
    int _size;
    T *pinned_ptr;

  public:
    pinnedVector() {}
    pinnedVector(int);
    pinnedVector(int, T &);
    ~pinnedVector();
    void alloc(int);
    void alloc(int, T &);
    int size();
    T &operator[](int);
    void copyToDevice(cudaVector<T> &, cudaStream_t &);
    void copyFromDevice(cudaVector<T> &, cudaStream_t &);
    void operator=(std::vector<T> &);
    void operator=(cudaVector<T> &);
    void copyTo(std::vector<T> &);
    T *getPointer();
};

template <class T> class pinnedScalar {

    T *pinned_ptr;

  public:
    pinnedScalar(T &val) {
        cudaHostAlloc((void **)&pinned_ptr, sizeof(T), cudaHostAllocDefault);
        *pinned_ptr = val;
    }

    pinnedScalar() {
        cudaHostAlloc((void **)&pinned_ptr, sizeof(T), cudaHostAllocDefault);
    }

    void copyToDevice(cudaScalar<T> &dev, cudaStream_t &stream) {
        cudaMemcpyAsync(dev.getPointer(), pinned_ptr, sizeof(T),
                        cudaMemcpyHostToDevice, stream);
    }

    void copyFromDevice(cudaScalar<T> &dev, cudaStream_t &stream) {
        cudaMemcpyAsync(pinned_ptr, dev.getPointer(), sizeof(T),
                        cudaMemcpyDeviceToHost, stream);
    }

    void operator=(T v) { *pinned_ptr = v; }

    void copyTo(T &v) { v = *pinned_ptr; }

    T *getPointer() { return pinned_ptr; }
};

#endif
