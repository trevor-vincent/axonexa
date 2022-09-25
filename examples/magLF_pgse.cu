#define WARP_SIZE 32
#define USE_DOUBLE
//#define SPECULAR_REFLECTION
//#define USE_RELAXATION
//#define GAMMA 267500.0 // ms^-1 * T^-1
#define GAMMA (-203789.0)
#define PI 3.1415926535897932384626433832795
//#define USE_RELAXATION

#include <assert.h>
#include <cmath>
#include <cuda.h>
#include <curand_kernel.h>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <vector>
// #define DORP_NORELAXATION
// #define DORP_NORELAXATION

//#define D_NOT_ZERO
#if defined USE_DOUBLE
typedef double real;
//#define EPSILON 1e-14 //good for small radii
#define EPSILON 5e-12 // good for large radii

#else
typedef float real;
#define EPSILON 1e-6

#endif

using namespace std;

#include "CPUkernels.cuh"
#include "Sphere.cuh"
#include "bfunctors.cuh"
#include "blochdiff.cuh"
#include "boundaryCheck.cuh"
#include "compare.cuh"
#include "cudaVector.cu"
#include "cudaVector.cuh"
#include "cylinderXY.cuh"
#include "empty.cuh"
#include "gfunctors.cuh"
#include "kernelDEBUG.cuh"
#include "kernelLattice.cuh"
#include "kernelMag.cuh"
#include "kernelPhase.cuh"
#include "kernelSetup.cuh"
#include "kernelWC.cuh"
#include "lattice.cuh"
#include "magAcquisition.cuh"
#include "magAcquisitionStream.cuh"
#include "misc.cuh"
#include "phaseAcquisition.cuh"
#include "phaseAcquisitionStream.cuh"
#include "pinnedVector.cu"
#include "pinnedVector.cuh"
#include "plane.cuh"
#include "simuparams.cuh"
#include "substrate.cuh"
#include "timer.cuh"
#include "vector3.cuh"

int main() {

    cudaFuncSetCacheConfig("updateWalkersMag", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("setup_kernel", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("_functionReduceAtom", cudaFuncCachePreferShared);

    int number_of_particles = 57344; // needs to be a factor of two
    real D = 0.0;
    real timestep = .000001;

    real speed = sqrt(6.0 * D / timestep);
    real rad = 24.2;
    // real rad = .005;
    int threads = 128;
    int blocks = number_of_particles / threads;

    int numOfAcq = 100;
    real B0 = .001;

    magAcquisitionStream<QSC_PGSE> pas(number_of_particles);
    // magAcquisitionStream<QSC_GRE> pas(number_of_particles);
    real gradient_duration = 1.25; //= echo_time*.25;
    real gradient_spacing = 0.2;   // echo_time - 2.0*gradient_duration;
    std::cout << "gradient_duration = " << gradient_duration << std::endl;
    std::cout << "gradient_spacing = " << gradient_spacing << std::endl;

    for (int i = 0; i < numOfAcq; i++) {

        // real echo_time = 2.5;
        real echo_time = 2.0 * gradient_duration + gradient_spacing;
        QSC_PGSE qsc(B0, i * .000006 / ((real)numOfAcq), gradient_duration,
                     gradient_spacing, timestep, .000001, timestep, echo_time);
        // QSC_GRE
        // qsc(B0,i*.000006/((real)numOfAcq),gradient_duration,gradient_spacing);
        echo_time += 2.0 * timestep;
        int number_of_timesteps = (int)(echo_time / timestep);
        // magAcquisition<QSC_GRE> pa(1,number_of_timesteps,i*time(NULL));
        magAcquisition<QSC_PGSE> pa(1, number_of_timesteps, i * time(NULL));
        pa.addMeasurement(qsc);
        pas.addAcquisition(pa);
    }

    // pas.runCPUAcquisition(0, rad, timestep,  c, speed);
    Cylinder_XY cylinder(0.0, 0.0, rad, 0.0, 0.0, D, 0, 0.0);
    Vector3 initialM(1.0, 0.0, 0.0);

    std::vector<int> plan(3);
    plan[0] = 0;
    plan[1] = numOfAcq;
    plan[2] = numOfAcq;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    pas.runAcquisitionStream(cylinder, initialM, timestep, blocks, threads, 1,
                             plan, numOfSMPerDevice);
    // for (int i = 0; i < numOfAcq; i++){
    // pas.runCPUAcquisition(i,cylinder, initialM,  timestep);
    // }

    for (int i = 0; i < numOfAcq; i++) {
        cout << setprecision(20) << endl;
        cout << pas.getAcquisition(i).getFieldFunctors()[0].G;
        cout << " ";
        cout << pas.getAcquisition(i).getFieldFunctors()[0].G *
                    pas.getAcquisition(i).getFieldFunctors()[0].G;
        cout << " ";
        cout << (Vector3(pas.getAcquisition(i).getMx()[0],
                         pas.getAcquisition(i).getMy()[0], 0.0) *
                 (1.0 / number_of_particles))
                    .magnitude();
    }
}
