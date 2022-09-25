#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
//#define USE_RELAXATION

#include <algorithm>
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
#include "RPSinitializer.h"
#include "Sphere.cuh"
#include "bfunctors.cuh"
#include "blochdiff.cuh"
#include "boundaryCheck.cuh"
#include "compare.cuh"
#include "cudaVector.cu"
#include "cudaVector.cuh"
#include "cylinderXY.cuh"
#include "deviates.h"
#include "empty.cuh"
#include "gamma.h"
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
#include "nr3.h"
#include "phaseAcquisition.cuh"
#include "phaseAcquisitionStream.cuh"
#include "pinnedVector.cu"
#include "pinnedVector.cuh"
#include "plane.cuh"
#include "ran.h"
#include "simuparams.cuh"
#include "substrate.cuh"
#include "timer.cuh"
#include "vector3.cuh"

int main() {

    cudaFuncSetCacheConfig("updateWalkersMag", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("setup_kernel", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("_functionReduceAtom", cudaFuncCachePreferShared);

    int number_of_particles = 114688; // needs to be a factor of two
    real timestep = .001;

    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<RectGrad> pas(number_of_particles);

    int NOI = 6;
    int NOM = 20;

    real gradient_duration = 10;

    for (int j = 0; j < NOI; j++) {
        real gradient_spacing = gradient_duration + j * 10;
        real echo_time = 2.0 * gradient_duration + gradient_spacing;
        int number_of_timesteps = (int)(echo_time / timestep);
        phaseAcquisition<RectGrad> pa(NOM, number_of_timesteps,
                                      number_of_particles, time(NULL));

        for (int i = 0; i < NOM; i++) {
            real G = 8. * (real(i) / real(NOM)) * 20 * 0.0000025;
            RectGrad cosGRAD(G, gradient_duration, gradient_spacing);
            pa.addMeasurement(cosGRAD);
        }

        pas.addAcquisition(pa);
    }

    real D_extra = 2.5E-6;
    real D_intra = 1.0E-6;
    real T2_i = 200;
    real T2_e = 200;

    real cube_length = .0365;
    real alpha = .00691;     // alpha = k
    real beta = 1.0 / 2.331; // beta = 1/theta
    real radmin = 2.0 * sqrt(6.0 * D_extra * timestep);

    std::vector<Cylinder_XY> basis;
    Lattice lattice(cube_length, cube_length, cube_length, T2_e, 0.0, D_extra,
                    100);
    RPSLatticeInitializer<Cylinder_XY> rpsli(lattice, 0);
    rpsli.gammaRadialDist(624124, alpha, beta, radmin, cube_length / 7.5);
    rpsli.uniformCenterDist(1208993908, 500);
    rpsli.setRegions();
    rpsli.correctEdges();
    lattice =
        rpsli.lat; // needed to reinitialize basis size (since it was
                   // initialized to 100 and there will be > 100 cylinders).

    // if (rpsli.basis.size() != 100){std::cout << " Basis Size Does not equal
    // 100 " << std::endl;} std::cout << "lattice basis size = " <<
    // lattice.getBasisSize() << std::endl;
    for (int i = 0; i < rpsli.basis.size(); i++) {
        basis.push_back(
            Cylinder_XY(0.0, 0.0, 0.0, T2_i, 0.0, D_intra, i + 1, 0.0));
        basis[i].setRadius(rpsli.basis[i].getRadius());
        basis[i].setCenter(rpsli.basis[i].getCenter());
        basis[i].setEPS((1E-13));
        basis[i].setRegion(rpsli.basis[i].getRegion());
        std::cout << rpsli.basis[i].getCenter() << "  "
                  << rpsli.basis[i].getRadius() << " " << std::endl;
    }

    std::cout << " Packing Fraction = " << rpsli.getVI() << std::endl;

    std::vector<int> plan(3);
    plan[0] = 0;
    plan[1] = NOI;
    plan[2] = NOI;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    pas.runAcquisitionStreamLattice(&basis[0], lattice, timestep, blocks,
                                    threads, 1, plan, numOfSMPerDevice);

    std::cout << std::endl << " Signals " << std::endl;

    for (int i = 0; i < NOI; i++) {
        for (int j = 0; j < NOM; j++) {
            std::cout << setprecision(20);
            std::cout << pas.getAcquisition(i).getGradientFunctors()[j].getTau()
                      << " ";
            std::cout << pas.getAcquisition(i).getGradientFunctors()[j].getG()
                      << " ";
            std::cout << pas.getAcquisition(i).getMx()[j] << " "
                      << pas.getAcquisition(i).getMy()[j] << " ";
            std::cout << std::endl;
        }
    }
}
