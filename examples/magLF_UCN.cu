#define WARP_SIZE 32
#define USE_DOUBLE
//#define SPECULAR_REFLECTION
//#define USE_RELAXATION
//#define GAMMA 267500.0 // ms^-1 * T^-1
// define GAMMA (-203789.0)
#define GAMMA 73997.0
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
#define ROTATE_MULTISTEP

// #define D_NOT_ZERO
#if defined USE_DOUBLE
typedef double real;
//#define EPSILON 1e-14  //good for small radii
//#define EPSILON 5e-12  //good for D=20, r ~ 25
#define EPSILON 5e-10 // good for D=880, r ~ 150
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
    real D = 880;
    real timestep = .001;

    real rad = 150;
    int threads = 128;
    int blocks = number_of_particles / threads;

    int numOfAcq = 1000;
    real B0 = .000001;
    std::vector<real> gradients(10);

    gradients[0] = 1E-8;
    gradients[1] = 5E-9;
    gradients[2] = 1E-9;
    gradients[3] = 5E-10;
    gradients[4] = 1E-10;
    gradients[5] = 5E-11;
    gradients[6] = 1E-11;
    gradients[7] = 5E-12;
    gradients[8] = 1E-12;
    gradients[9] = 5E-13;

    magAcquisitionStream<UCN_lowfield> pas(number_of_particles);

    for (int i = 0; i < numOfAcq; i++) {

        // real echo_time = 2.5;
        int number_of_timesteps = (i + 1) * 100;
        magAcquisition<UCN_lowfield> pa(gradients.size(), number_of_timesteps,
                                        i * time(NULL));

        for (int j = 0; j < gradients.size(); j++) {
            UCN_lowfield qsc(B0, gradients[j], timestep, real(.000001),
                             timestep);
            pa.addMeasurement(qsc);
        }

        pas.addAcquisition(pa);
    }

    // pas.runCPUAcquisition(0, rad, timestep,  c, speed);
    Cylinder_XY cylinder(0.0, 0.0, rad, 0.0, 0.0, D, 0, 0.0, EPSILON);
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
    cout << setprecision(20) << endl;

    for (int j = 0; j < gradients.size(); j++) {
        for (int i = 0; i < numOfAcq; i++) {
            std::cout << std::endl;
            cout << pas.getAcquisition(i).getFieldFunctors()[j].G;
            cout << " ";
            cout << (i + 1) * 100 * timestep;
            cout << " ";
            cout << pas.getAcquisition(i).getMx()[j];
            cout << " ";
            cout << pas.getAcquisition(i).getMy()[j];
            cout << " ";
            cout << pas.getAcquisition(i).getMz()[j];
        }
    }
}
