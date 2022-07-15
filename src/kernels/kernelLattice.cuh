template <class GradientType, class Basis>
__global__ void updateWalkersLatticePhase(

    const SimuParamsPhase *par, const Lattice *lat, const Basis *basis,
    const GradientType *G, curandState *globalState, real *T2factors,
    real *phase

)

{

    const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    curandState localState = globalState[tid];

    // accumTime is the time spent by a walker in a region before transmission.
    // accumTime is equal to the timestep if the particle never gets
    // transmitted, but it is in general smaller than the timestep. It is needed
    // for calculating T2, which requires the time spent in each region during a
    // timestep.

    real phi, theta;
#if defined USE_RELAXATION
    real T2factor = 0.0;
#endif

    Vector3 r, ri;
    real speed;
#if defined USE_INITIALIZE_IN_REGION
#if USE_INITIALIZE_IN_REGION == 0
    lat->initializeInRegion(basis, localState, r, 0);
#else
    basis[USE_INITIALIZE_IN_REGION - 1].randUnif(r, localState);
#endif
#else
    lat->initializeUniformly(r, localState);
#endif
    Vector3 r_unbounded = r;

    for (int i = 0; i < par->measurements; i++) {
        // phase[tid + i*par->number_of_particles] = par->phase_initial;
        phase[tid + i * par->number_of_particles] =
            GAMMA * (G[i](0.0) * r) * par->timestep;
    }

    for (int i = 1; i < par->steps; i++) {

        ri = r;

        speed = sqrt(6.0 * lat->getD(basis, r) / par->timestep);

        phi = 2.0 * PI * curand_uniform(&localState);
        theta = acos(2.0 * curand_uniform(&localState) - 1);

        r += Vector3(speed * par->timestep * sin(theta) * cos(phi),
                     speed * par->timestep * sin(theta) * sin(phi),
                     speed * par->timestep * cos(theta));

// Note that the permeable case also includes the impermeable one (just set P =
// 0 and set PERMEABLE). However, to lower register usage, we have opted to have
// a different set of functions for the impermeable case as well. This allows
// for speedier calculations when all geometry elements have impermeable
// boundaries.
#if defined SPECULAR_REFLECTION
#if defined USE_PERMEABLE
        real accumtime = 0.0;
        boundaryCheckPermeable(lat, basis, ri, r, accumtime, par->timestep,
                               localState);
#else
        boundaryCheckImpermeable(lat, basis, ri, r, speed, par->timestep);
#endif

#else // SIMPLE REJECTION

        if (lat->inRegion(basis, ri) != lat->inRegion(basis, r)) {
            r = ri;
        }

#endif

        r_unbounded += (r - ri);
        for (int j = 0; j < par->measurements; j++) {
            phase[tid + j * par->number_of_particles] +=
                GAMMA * (G[j](par->timestep * i) * r_unbounded) * par->timestep;
        }

#if defined USE_RELAXATION
#if defined USE_PERMEABLE
        real T2_i = lat->getT2(basis, ri);
        real T2_f = lat->getT2(basis, r);
        T2factor += (accumtime / T2_i) + ((par->timestep - accumtime) / T2_f);
#else
        real T2 = lat->getT2(basis, ri);
        T2factor += par->timestep / T2;
#endif
#endif

        lat->correctBoundary(r);
    }

#if defined USE_RELAXATION
    T2factors[tid] = -T2factor;
#endif
}
