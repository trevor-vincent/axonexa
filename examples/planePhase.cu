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

    int number_of_particles = std::pow(2,12); // needs to be a factor of two
    real timestep = .001;

    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<CosGFunc> pas(number_of_particles);

    int NOI = 10;
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

    std::vector<real> lengths(4);
    lengths[0] = .001;
    lengths[1] = .002;
    lengths[2] = .005;
    lengths[3] = .010;

    for (int r = 0; r < 4; r++) {

        real D_extra = 2.5E-6;
        real D_intra = 1.0E-6;
        real T2_i = 200;
        real T2_e = 200;

        real l = lengths[r];
        std::cout << " PLANE LENGTH = " << l << std::endl;
        Plane cyls[2];
        Lattice lattice(l, l, l, T2_e, 0.0, D_extra, 2);
        cyls[0] = Plane(Vector3(0.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0));
        cyls[1] = Plane(Vector3(l, 0.0, 0.0), Vector3(-1.0, 0.0, 0.0));

        // Empty nothing[1];
        // Lattice lattice(d, 2.0*d*0.86602540378443864676372317075294, d, T2_e,
        // 0.0, D_extra); nothing[0] = Empty();

        // for (int i = 0; i < NOI; i++){
        // pas.runAcquisitionLattice(i,cyls, lattice, timestep, blocks, threads,
        // 14);
        // }

        std::vector<int> plan(3);
        plan[0] = 0;
        plan[1] = NOI;
        plan[2] = NOI;
        std::vector<int> numOfSMPerDevice(1);
        numOfSMPerDevice[0] = 1;
        numOfSMPerDevice[1] = 2;

        for (int i = 0; i < 3; i++) {
            pas.runAcquisitionStreamLattice(cyls, lattice, timestep, blocks,
                                            threads, 1, plan, numOfSMPerDevice);
            pas.calcEveryADC();
            pas.changeSeeds(time(NULL) * (i + 1));
        }

        for (int i = 0; i < NOI; i++) {
            std::cout << setprecision(20);
            std::cout
                << pas.getAcquisition(i).getGradientFunctors()[0].getFreq()
                << " ";
            std::cout << pas.getAcquisition(i).getADC() << std::endl;
        }
        pas.flushADC();
    }
}
