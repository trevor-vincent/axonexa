// accumulates the b value

template <class SEQUENCE> double bvalue(SEQUENCE *seq, double dt) {

    // cout << "b value CALCULATION " << endl;
    int num_of_steps = (int)ceil(seq->get_TE() / dt);
    // cout << "num of steps = " << num_of_steps << endl;

    double b = 0.0, Gint;

    for (int i = 0; i < num_of_steps; i++) {

        Gint = gradint(i * dt, seq, dt);
        b += Gint * Gint * dt;
    }

    return b;
}

template <class SEQUENCE> double gradint(double t, SEQUENCE *seq, double dt) {

    int num_of_steps = (int)ceil(t / dt);
    double Gint = 0.0;

    for (int i = 0; i < num_of_steps; i++) {
        Gint += GAMMA * (seq->gradient(i * dt)) * dt;
    }

    return Gint;
}

// void getGPUready(){

// }

// __global__ void bvalue_SIN(float* binc, int sizeofbinc, float grad_strength,
// float grad_duration, float grad_spacing ){

// int tid = blockIdx.x*blockDim.x + threadIdx.x;
// float GAMMA = 267522.0

// while (tid_x < sizeofbinc){

// for (int i = 0; i < tid_x; i++){

// binc[tid] += SIN_OGSE_GRAD(i*dt, period, grad_strength, grad_duration,
// grad_spacing)*GAMMA*timestep;

// }

// binc[tid] *= binc[tid_x];
// binc[tid] *= timestep;
// tid += blockDim.x * gridDim.x;

// }

// tid = blockIdx.x*blockDim.x + threadIdx.x;

// if (sizeofbinc % 2 == 0){
// int index = sizeofbinc/2;

// while (tid < index){

// b[tid] += b[tid + index];
// index/=2;

// }

// }

// else {

// }

// }

// __device__ float SIN_OGSE_GRAD(float time, float period, float grad_strength,
// float grad_duration, float grad_spacing){

// float Pie = 3.141592654f;

// if (time < grad_duration){

// return grad_strength*sin(2*pie*(time)/period);

// }

// else if ( (fabs(time - grad_spacing - grad_duration) < timestep/2.0f || time
// > grad_spacing + grad_duration) && time < 2*grad_spacing + grad_duration){

// return -grad_strength*sin(2*pie*(time - grad_spacing -
// grad_duration)/period);
// }

// else {
// return 0;
// }

// }
