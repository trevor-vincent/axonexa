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

    phaseAcquisitionStream<CosGFunc> pas(number_of_particles);
    phaseAcquisitionStream<CosGFunc> pas1;
    phaseAcquisitionStream<CosGFunc> pas2;
    phaseAcquisitionStream<CosGFunc> pas3;
    phaseAcquisitionStream<CosGFunc> pas4;
    phaseAcquisitionStream<CosGFunc> pas5;
    phaseAcquisitionStream<CosGFunc> pas6;
    phaseAcquisitionStream<CosGFunc> pas7;
    phaseAcquisitionStream<CosGFunc> pas8;
    phaseAcquisitionStream<CosGFunc> pas9;
    phaseAcquisitionStream<CosGFunc> pas10;

    int NOI = 20;
    int NOM = 20;

    real gradient_duration = 20;
    real gradient_spacing = 2.0;
    real echo_time = 2.0 * gradient_duration + gradient_spacing;
    int number_of_timesteps = (int)(echo_time / timestep);
    phaseAcquisition<CosGFunc> pa(NOM * NOI, number_of_timesteps,
                                  number_of_particles, time(NULL));

    for (int j = 0; j < NOI; j++) {
        for (int i = 0; i < NOM; i++) {
            int N = 1 + 10 * j;
            real G = 8. * i * 0.0000025 * N;
            CosGFunc cosGRAD(G, gradient_duration, gradient_spacing, N,
                             Vector3(1.0, 0.0, 0.0));
            pa.addMeasurement(cosGRAD);
        }
    }

    pas.addAcquisition(pa);

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
    plan[1] = NOI * NOM;
    plan[2] = NOI * NOM;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    pas1 = pas;
    pas1.getAcquisition(0).getSeed() *= 2;
    pas2 = pas;
    pas2.getAcquisition(0).getSeed() *= 3;
    pas3 = pas;
    pas3.getAcquisition(0).getSeed() *= 4;
    pas4 = pas;
    pas4.getAcquisition(0).getSeed() *= 5;
    pas5 = pas;
    pas4.getAcquisition(0).getSeed() *= 6;
    pas6 = pas;
    pas4.getAcquisition(0).getSeed() *= 7;
    pas7 = pas;
    pas4.getAcquisition(0).getSeed() *= 8;
    pas8 = pas;
    pas4.getAcquisition(0).getSeed() *= 9;
    pas9 = pas;
    pas4.getAcquisition(0).getSeed() *= 10;

    pas.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                              14);
    pas1.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas2.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas3.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas4.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas5.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas6.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas7.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas8.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas9.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);

    std::cout << std::endl << " Signals " << std::endl;

    for (int j = 0; j < NOI * NOM; j++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(0).getGradientFunctors()[j].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(0).getGradientFunctors()[j].getG()
                  << " ";
        std::cout << pas.getAcquisition(0).getMx()[j] << " "
                  << pas.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas1.getAcquisition(0).getMx()[j] << " "
                  << pas1.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas2.getAcquisition(0).getMx()[j] << " "
                  << pas2.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas3.getAcquisition(0).getMx()[j] << " "
                  << pas3.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas4.getAcquisition(0).getMx()[j] << " "
                  << pas4.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas5.getAcquisition(0).getMx()[j] << " "
                  << pas5.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas6.getAcquisition(0).getMx()[j] << " "
                  << pas6.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas7.getAcquisition(0).getMx()[j] << " "
                  << pas7.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas8.getAcquisition(0).getMx()[j] << " "
                  << pas8.getAcquisition(0).getMy()[j] << " ";
        std::cout << pas9.getAcquisition(0).getMx()[j] << " "
                  << pas9.getAcquisition(0).getMy()[j] << " ";
        std::cout << std::endl;
    }
}
