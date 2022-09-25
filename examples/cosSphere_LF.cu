#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
//#define USE_RELAXATION
#define D_NOT_ZERO

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
    real B0 = .00001;

    for (int r = 0; r < 5; r++) {

        B0 *= (r > 0) * 10.0 + (r == 0) * 1.0;
        std::cout << " B0 = " << B0 << std::endl;

        magAcquisitionStream<QSC_COS_GRE> pas(number_of_particles);
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
            real gradient_spacing = .5;
            real echo_time = 2.0 * gradient_duration + gradient_spacing;
            int number_of_timesteps = (int)(echo_time / timestep);
            magAcquisition<QSC_COS_GRE> pa(NOM, number_of_timesteps,
                                           j * time(NULL));

            for (int i = 0; i < NOM; i++) {

                int N = 1 + j;
                real G = i * 0.0000025 * N;

                QSC_COS_GRE qscCOS(B0, G, gradient_duration, gradient_spacing,
                                   N);
                pa.addMeasurement(qscCOS);
            }
            pas.addAcquisition(pa);
        }

        real rad3 = .0025;
        Sphere spher3(0.0, 0.0, 0.0, rad3, 0.0, 0.0, D, 0, 0.0);

        std::vector<int> plan(3);
        plan[0] = 0;
        plan[1] = NOI;
        plan[2] = NOI;
        std::vector<int> numOfSMPerDevice(1);
        numOfSMPerDevice[0] = 14;
        numOfSMPerDevice[1] = 2;
        Vector3 initialM(1.0, 0.0, 0.0);

        CPUtimer ctime;
        ctime.start();
        for (int i = 0; i < 3; i++) {
            pas.runAcquisitionStream(spher3, initialM, timestep, blocks,
                                     threads, 1, plan, numOfSMPerDevice);
            pas.calcEveryADC();
            pas.changeSeeds(time(NULL) * (i + 1));
        }
        ctime.stop();
        for (int i = 0; i < NOI; i++) {
            std::cout << setprecision(20);
            std::cout << pas.getAcquisition(i).getFieldFunctors()[0].getFreq()
                      << " ";
            std::cout << pas.getAcquisition(i).getADC() << " ";
            std::cout << pas.getAcquisition(i).getADCError() << std::endl;
        }
        ctime.display();
    }
}
