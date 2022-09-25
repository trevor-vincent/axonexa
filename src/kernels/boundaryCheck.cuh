#include "lattice.cuh"

template <class Basis>
__device__ void boundaryNormal(Vector3 &ri, Vector3 &r, const real currentSpeed,
                               const Basis *basis, const real timestep) {

    real v = 0.0;
    real accumtime = 0.0;

    while (basis->intersect(ri, r, v)) {

        r = (r - ri) * v + ri;

        // determine incident vector
        Vector3 I = (r - ri);
        real Imag = I.magnitude();

        // calculate accumulated time
        accumtime += Imag / currentSpeed;

        // determine normal for reflection calculation
        Vector3 n = basis->getNormal(r);

        // set last position to boundary position before rflection
        ri = r;

        // reflect, reflection travel time is the remaining time (timestep -
        // accumtime)
        r += (I - n * 2.0 * (n * I)) * (currentSpeed / Imag) *
             (timestep - accumtime);
    }
}

template <class Basis>
__device__ void boundaryCheckImpermeable(const Lattice *lat, const Basis *basis,
                                         const Vector3 &ri, Vector3 &r,
                                         const real currentSpeed,
                                         const real timestep) {

    real v = 0.0;
    real accumtime = 0.0;
    Vector3 n;
    Vector3 riCurrent = ri;

    while (lat->intersectionCheckImpermeable(basis, riCurrent, r, v, n)) {

        r = (r - riCurrent) * v + riCurrent;

        // determine incident vector
        Vector3 I = (r - riCurrent);
        real Imag = I.magnitude();

        // calculate accumulated time
        accumtime += Imag / currentSpeed;

        // set last position to boundary position before rflection
        riCurrent = r;

        // reflect, reflection travel time is the remaining time (timestep -
        // accumtime)
        r += (I - n * 2.0 * (n * I)) * (currentSpeed / Imag) *
             (timestep - accumtime);
    }
}

template <class Basis>
__device__ void boundaryCheckPermeable(const Lattice *lat, const Basis *basis,
                                       const Vector3 &ri, Vector3 &r,
                                       real accumtime, const real timestep,
                                       curandState &localState) {

    real v = 0.0;
    Vector3 n;
    bool transmitted = false;
    real permeability;
    real currentSpeed = sqrt(6.0 * lat->getD(basis, ri) / timestep);
    Vector3 riCurrent = ri;

    while (lat->intersectionCheckPermeable(basis, riCurrent, r, v, n,
                                           permeability) &&
           transmitted == false) {

        real transmissionSpeed = sqrt(6.0 * lat->getD(basis, r) / timestep);

        r = (r - riCurrent) * v + riCurrent;

        // determine incident vector
        Vector3 I = (r - riCurrent);
        real Imag = I.magnitude();

        // calculate accumulated time
        accumtime += Imag / currentSpeed;

        // reflect
        if (curand_uniform(&localState) > permeability * 4.0 / currentSpeed) {

            // set last position to boundary position before rflection
            riCurrent = r;

            // reflect, reflection travel time is the remaining time (timestep -
            // accumtime)
            r += (I - n * 2.0 * (n * I)) * (currentSpeed / Imag) *
                 (timestep - accumtime);

        }

        else {

            r += (I / Imag) * transmissionSpeed * (timestep - accumtime);
            transmitted = true;
        }
    }
}
