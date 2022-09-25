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
#include "simuparams.cuh"
#include "substrate.cuh"
#include "timer.cuh"
#include "vector3.cuh"

int main() {

    cudaFuncSetCacheConfig("updateWalkersMag", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("setup_kernel", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("_functionReduceAtom", cudaFuncCachePreferShared);

    int number_of_particles = 114688; // needs to be a factor of two
    real D = 2.5E-6;
    // real D = 1E-6;
    real timestep = .001;
    real speed = sqrt(6.0 * D / timestep);

    // real rad = .005;
    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<CosGFunc> pas(number_of_particles);
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

    int NOI = 100;
    int NOM = 2;

    for (int j = 0; j < NOI; j++) {

        real gradient_duration = 10;
        real gradient_spacing = 2.0;
        real echo_time = 2.0 * gradient_duration + gradient_spacing;
        int number_of_timesteps = (int)(echo_time / timestep);
        phaseAcquisition<CosGFunc> pa(NOM, number_of_timesteps,
                                      number_of_particles, j * time(NULL));

        for (int i = 0; i < NOM; i++) {

            int N = 1 + j;
            real G = i * 0.0000025 * N;

            CosGFunc cosGRAD(G, gradient_duration, gradient_spacing, N,
                             Vector3(1.0, 0.0, 0.0));
            pa.addMeasurement(cosGRAD);
        }
        pas.addAcquisition(pa);
    }

    real rad1 = .0005;
    real rad2 = .001;
    real rad3 = .0025;
    real rad4 = .005;

    Cylinder_XY spher1(0.0, 0.0, rad1, 0.0, 0.0, D, 0, 0.0);
    Cylinder_XY spher2(0.0, 0.0, rad2, 0.0, 0.0, D, 0, 0.0);
    Cylinder_XY spher3(0.0, 0.0, rad3, 0.0, 0.0, D, 0, 0.0);
    Cylinder_XY spher4(0.0, 0.0, rad4, 0.0, 0.0, D, 0, 0.0);

    std::vector<int> plan(3);
    plan[0] = 0;
    plan[1] = NOI;
    plan[2] = NOI;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    for (int i = 0; i < 3; i++) {
        pas.runAcquisitionStream(spher1, timestep, blocks, threads, 1, plan,
                                 numOfSMPerDevice);
        pas.calcEveryADC();
        pas.changeSeeds(time(NULL) * (i + 1));
    }
    std::cout << std::endl << std::endl;

    for (int i = 0; i < NOI; i++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(i).getADC() << std::endl;
    }
    pas.flushADC();

    std::cout << " NEXT Cylinder_XY " << std::endl;

    for (int i = 0; i < 3; i++) {
        pas.runAcquisitionStream(spher2, timestep, blocks, threads, 1, plan,
                                 numOfSMPerDevice);
        pas.calcEveryADC();
        pas.changeSeeds(time(NULL) * (i + 1));
    }
    std::cout << std::endl << std::endl;

    for (int i = 0; i < NOI; i++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(i).getADC() << std::endl;
    }
    pas.flushADC();

    std::cout << " NEXT Cylinder_XY " << std::endl;

    for (int i = 0; i < 3; i++) {
        pas.runAcquisitionStream(spher3, timestep, blocks, threads, 1, plan,
                                 numOfSMPerDevice);
        pas.calcEveryADC();
        pas.changeSeeds(time(NULL) * (i + 1));
    }
    std::cout << std::endl << std::endl;

    for (int i = 0; i < NOI; i++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(i).getADC() << std::endl;
    }
    pas.flushADC();

    std::cout << " NEXT Cylinder_XY " << std::endl;

    for (int i = 0; i < 3; i++) {
        pas.runAcquisitionStream(spher4, timestep, blocks, threads, 1, plan,
                                 numOfSMPerDevice);
        pas.calcEveryADC();
        pas.changeSeeds(time(NULL) * (i + 1));
    }

    std::cout << std::endl << std::endl;

    for (int i = 0; i < NOI; i++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(i).getADC() << std::endl;
    }
    pas.flushADC();
}
