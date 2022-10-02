#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
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

#if defined USE_DOUBLE
typedef double real;
#define EPSILON 1e-14

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
#include "simuParams.cuh"
#include "substrate.cuh"
#include "timer.cuh"
#include "vector3.cuh"

int main() {

    cudaFuncSetCacheConfig("updateWalkersMag", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("setup_kernel", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("_functionReduceAtom", cudaFuncCachePreferShared);

    int number_of_particles = 16384; // needs to be a factor of two
    real D = 20.0;
    // real D = 1E-6;
    real timestep = .001;

    real speed = sqrt(6.0 * D / timestep);
    real rad = 50;
    // real rad = .005;
    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<RectGrad> pas(number_of_particles);
    // int echo_time = 1;
    // int number_of_timesteps = (int) (echo_time/timestep);
    // real gradient_duration = echo_time*.25;
    // real gradient_spacing = echo_time - 2.0*gradient_duration;
    // QSC_GRE
    // qsc(.0006, 1.5811388300841896659994467722164e-4,gradient_duration,gradient_spacing);
    // magAcquisition<QSC_GRE>
    // pa(1,number_of_timesteps,number_of_particles,1234);
    // pa.addMeasurement(qsc);
    // pas.addAcquisition(pa);

    int numOfAcq = 20;

    for (int i = 0; i < numOfAcq; i++) {

        real echo_time = 1;
        int number_of_timesteps = (int)(echo_time / timestep);
        real gradient_duration = echo_time * .25;
        real gradient_spacing = echo_time - 2.0 * gradient_duration;
        RectGrad qsc(i * 1000000 / ((real)numOfAcq), gradient_duration,
                     gradient_spacing, Vector3(1.0, 0.0, 0.0));
        // std::cout << "echo_time = " << echo_time << std::endl;
        // std::cout << "gradient_duration = " << gradient_duration <<
        // std::endl; std::cout << "gradient_spacing = " << gradient_spacing <<
        // std::endl;

        // RECT_GRE qsc(0.0,
        // Vector3(i*7.0710678118654752440084436210485e-6/1000.0, 0.0,
        // 0.0),gradient_duration,gradient_spacing);
        phaseAcquisition<RectGrad> pa(1, number_of_timesteps,
                                      number_of_particles, i * time(NULL));
        pa.addMeasurement(qsc);
        pas.addAcquisition(pa);
    }

    // pas.runCPUAcquisition(0, rad, timestep,  c, speed);
    Cylinder_XY cylinder(0.0, 0.0, rad, 0.0, 0.0, D, 0);
    Vector3 initialM(1.0, 0.0, 0.0);
    Substrate sub(90, 90, .001);

    pas.runAcquisition(1, cylinder, sub, initialM, timestep, blocks, threads,
                       14);
    cout << pas.getAcquisition(1).getMx()[0] << endl;
    cout << pas.getAcquisition(1).getMy()[0] << endl;

    pas.runAcquisitionWC(1, cylinder, sub, initialM, timestep, blocks, threads);
    cout << pas.getAcquisition(1).getMx()[0] << endl;
    cout << pas.getAcquisition(1).getMy()[0] << endl;

    std::vector<int> plan(3);
    plan[0] = 0;
    plan[1] = numOfAcq;
    plan[2] = numOfAcq;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    pas.runAcquisitionStream(cylinder, sub, initialM, timestep, blocks, threads,
                             1, plan, numOfSMPerDevice);

    for (int i = 0; i < numOfAcq; i++) {
        cout << setprecision(20) << endl;
        cout << Vector3(pas.getAcquisition(i).getMx()[0],
                        pas.getAcquisition(i).getMy()[0], 0.0) *
                    (1.0 / number_of_particles);
    }
}
