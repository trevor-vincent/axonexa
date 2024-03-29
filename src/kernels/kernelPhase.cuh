#pragma once
#include "vector3.cuh"

template <class GradientType, class Basis>
__global__ void updateWalkersPhase(const SimuParamsPhase *par,
                                   const Basis *basis, const GradientType *G,
                                   curandState *globalState, real *phase)

{

    const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curandState localState = globalState[tid];

    real phi, theta;
    real speed = sqrt(6.0 * basis->getD() / par->timestep);
    real stepLength = speed * par->timestep;

    Vector3 r = basis[0].unifRand(localState);

    for (int i = 0; i < par->measurements; i++) {
        // phase[tid + i*par->number_of_particles] = par->phase_initial;
        phase[tid + i * par->number_of_particles] =
            GAMMA * (G[i](0.0) * r) * par->timestep;
    }

    for (int i = 1; i < par->steps; i++) {

        Vector3 ri = r;

#if defined LATTICE_PICKING

        r += Vector3(stepLength * randomSign(localState),
                     stepLength * randomSign(localState),
                     stepLength * randomSign(localState));

#else

        phi = 2.0 * PI * curand_uniform(&localState);
        theta = acos(2.0 * curand_uniform(&localState) - 1);

        r += Vector3(stepLength * sin(theta) * cos(phi),
                     stepLength * sin(theta) * sin(phi),
                     stepLength * cos(theta));
#endif

#if defined SPECULAR_REFLECTION

        boundaryNormal(ri, r, speed, basis, par->timestep);

#else

        if (!basis[0].inside(r)) {
            r = ri;
        }

#endif

        for (int j = 0; j < par->measurements; j++) {
            phase[tid + j * par->number_of_particles] +=
                GAMMA * (G[j](par->timestep * i) * r) * par->timestep;
        }
    }

    // globalState[tid] = localState;
}
