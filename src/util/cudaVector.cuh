#ifndef CUDAVECTOR_H_
#define CUDAVECTOR_H_

template <class X> class pinnedVector;

template <class T> class cudaScalar {

    T _v;
    T *dev_ptr;
    int device;
    bool alloc;
    bool copied;

  public:
    cudaScalar(T val) {
        _v = val;
        safe_cuda(cudaMalloc((void **)&dev_ptr, sizeof(T)));
        alloc = true;
        copied = false;
        device = 0;
    }

    cudaScalar() {
        alloc = false;
        copied = false;
        device = 0;
    }

    ~cudaScalar() {
        if (alloc) {
            safe_cuda(cudaSetDevice(device));
            safe_cuda(cudaFree(dev_ptr));
        }
    }

    void setDevice(int _device) { device = _device; }

    void copyToDevice() {
        if (alloc == false) {
            safe_cuda(cudaMalloc((void **)&dev_ptr, sizeof(T)));
            alloc = true;
        }
        safe_cuda(cudaMemcpy(dev_ptr, &_v, sizeof(T), cudaMemcpyHostToDevice));
        copied = true;
    }

    void copyFromDevice() {
        if (alloc == true) {
            safe_cuda(
                cudaMemcpy(&_v, dev_ptr, sizeof(T), cudaMemcpyDeviceToHost));
        }
    }

    void operator=(T &v) { _v = v; }

    void malloc() {

        safe_cuda(cudaMalloc((void **)&dev_ptr, sizeof(T)));
        alloc = true;
    }

    T *getPointer() { return dev_ptr; }

    T getValue() { return _v; }
};

template <class T> class cudaVector {

  private:
    std::vector<T> _v;
    T *dev_ptr;
    int device;
    bool alloc;
    bool copied;

  public:
    cudaVector();
    cudaVector(int);
    cudaVector(int, T);
    ~cudaVector();
    void copyToDevice();
    void copyFromDevice();
    void add(T &);
    void setDevice(int);
    void malloc(int);
    int size();
    void resize(int);
    T &operator[](int);
    void operator=(pinnedVector<T> &);
    T *getPointer();
    void copyTo(std::vector<T> &);
    void operator=(std::vector<T> &);
    void sum(cudaScalar<T> &, int, int, cudaStream_t);
    void sum(cudaVector<T> &, int, int, int, int, cudaStream_t);
    template <class Transform>
    void transformAndSum(cudaScalar<T> &, int, int, cudaStream_t);
    template <class Transform>
    void transformAndSum(cudaVector<T> &, int, int, int, int, cudaStream_t);
    template <class Transform>
    void transformAndSumTwoVectors(cudaVector<T> &, cudaVector<T> &, int, int,
                                   int, int, cudaStream_t);
};

#endif
