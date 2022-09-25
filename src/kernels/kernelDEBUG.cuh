#include "simuParams.cuh"
#include "substrate.cuh"

template <class FieldType, class Basis>
__global__ void updateWalkersMagDEBUG(

    const SimuParams *par, const Substrate *sub, const Basis *basis,
    const FieldType *B, curandState *globalState,

    real *Mx, real *My, real *Mz,

    real *xi, real *yi, real *zi,

    real *x, real *y, real *z

)

{

    const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curandState localState = globalState[tid];

    real phi, theta;
    real speed = sqrt(6.0 * basis->getD() / par->timestep);

    // printf ( " timestep = %.6f \n" , par->timestep);
    // printf ( " measurements = %d \n" , par->measurements);
    // printf ( " number_of_particles = %d \n" , par->number_of_particles);
    // printf ( " steps = %d \n" , par->steps);
    // printf ( " c = %.6f \n" , sub->getC() );

    Vector3 r = basis[0].unifRand(localState, sub);
    xi[tid] = r.x;
    yi[tid] = r.y;
    zi[tid] = r.z;

    for (int i = 0; i < par->measurements; i++) {

        Mx[tid + i * par->number_of_particles] = par->mx_initial;
        My[tid + i * par->number_of_particles] = par->my_initial;
        Mz[tid + i * par->number_of_particles] = par->mz_initial;

        // printf ( " Kernel Minitial = %.6f %.6f %.6f \n",  Mx[tid +
        // i*par->number_of_particles], My[tid + i*par->number_of_particles],
        // Mz[tid + i*par->number_of_particles] );

        updateMag_noT2_rotate(&Mx[tid + i * par->number_of_particles],
                              &My[tid + i * par->number_of_particles],
                              &Mz[tid + i * par->number_of_particles],
                              B[i](r, (real)0.0), par->timestep);
    }

    for (int i = 1; i < par->steps; i++) {

        Vector3 ri = r;

        phi = 2.0 * PI * curand_uniform(&localState);
        theta = acos(2.0 * curand_uniform(&localState) - 1);

        r += Vector3(speed * par->timestep * sin(theta) * cos(phi),
                     speed * par->timestep * sin(theta) * sin(phi),
                     speed * par->timestep * cos(theta));

#if defined SPECULAR_REFLECTION

        boundaryNormal(ri, r, speed, basis, par->timestep);

#else

        if (!basis[0].inside(r)) {
            r = ri;
        }

#endif

        for (int j = 0; j < par->measurements; j++) {

            updateMag_noT2_rotate(&Mx[tid + j * par->number_of_particles],
                                  &My[tid + j * par->number_of_particles],
                                  &Mz[tid + j * par->number_of_particles],
                                  B[j](r, i * par->timestep), par->timestep);
            // printf ( " Kernel Mafter = %.6f %.6f %.6f \n",  Mx[tid +
            // j*par->number_of_particles], My[tid +
            // j*par->number_of_particles], Mz[tid + j*par->number_of_particles]
            // );
        }
    }

    // globalState[tid] = localState;

    x[tid] = r.x;
    y[tid] = r.y;
    z[tid] = r.z;
}
