#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
#define USE_INITIALIZE_IN_REGION 0
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

    int number_of_particles = 57344; // needs to be a factor of two
    real timestep = .001;

    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<CosGFunc> pas(number_of_particles);

    int NOI = 10;
    int NOM = 30;

    for (int j = 0; j < NOI; j++) {

        real gradient_duration = 20;
        real gradient_spacing = 2.0;
        real echo_time = 2.0 * gradient_duration + gradient_spacing;
        int number_of_timesteps = (int)(echo_time / timestep);
        phaseAcquisition<CosGFunc> pa(NOM, number_of_timesteps,
                                      number_of_particles, j * time(NULL));

        for (int i = 0; i < NOM; i++) {
            int N = 1 + 10 * j;
            real G = i * 0.0000025 * N;

            CosGFunc cosGRAD(G, gradient_duration, gradient_spacing, N,
                             Vector3(1.0, 0.0, 0.0));
            pa.addMeasurement(cosGRAD);
        }

        pas.addAcquisition(pa);
    }

    real radius = .0015;
    real D_extra = 2.5E-6;
    real D_intra = 1.0E-6;
    real f = .4;
    real a = pow((16. / 3.) * PI * radius * radius * radius / f, 1. / 3.);

    std::vector<Sphere> basis(14);
    Lattice lattice(a, a, a, 0.0, 0.0, D_extra, 14);
    basis[0] = Sphere(0.0 * a, 0.0 * a, 0.0 * a, radius, 0.0, 0.0, D_intra, 1,
                      0.0, EPSILON);
    basis[1] = Sphere(1.0 * a, 0.0 * a, 0.0 * a, radius, 0.0, 0.0, D_intra, 2,
                      0.0, EPSILON);
    basis[2] = Sphere(1.0 * a, 0.0 * a, 1.0 * a, radius, 0.0, 0.0, D_intra, 3,
                      0.0, EPSILON);
    basis[3] = Sphere(0.0 * a, 0.0 * a, 1.0 * a, radius, 0.0, 0.0, D_intra, 4,
                      0.0, EPSILON);
    basis[4] = Sphere(0.0 * a, 1.0 * a, 0.0 * a, radius, 0.0, 0.0, D_intra, 5,
                      0.0, EPSILON);
    basis[5] = Sphere(1.0 * a, 1.0 * a, 0.0 * a, radius, 0.0, 0.0, D_intra, 6,
                      0.0, EPSILON);
    basis[6] = Sphere(1.0 * a, 1.0 * a, 1.0 * a, radius, 0.0, 0.0, D_intra, 7,
                      0.0, EPSILON);
    basis[7] = Sphere(0.0 * a, 1.0 * a, 1.0 * a, radius, 0.0, 0.0, D_intra, 8,
                      0.0, EPSILON);
    basis[8] = Sphere(0.5 * a, 0.0 * a, 0.5 * a, radius, 0.0, 0.0, D_intra, 9,
                      0.0, EPSILON);
    basis[9] = Sphere(1.0 * a, 0.5 * a, 0.5 * a, radius, 0.0, 0.0, D_intra, 10,
                      0.0, EPSILON);
    basis[10] = Sphere(0.5 * a, 1.0 * a, 0.5 * a, radius, 0.0, 0.0, D_intra, 11,
                       0.0, EPSILON);
    basis[11] = Sphere(0.0 * a, 0.5 * a, 0.5 * a, radius, 0.0, 0.0, D_intra, 12,
                       0.0, EPSILON);
    basis[12] = Sphere(0.5 * a, 0.5 * a, 1.0 * a, radius, 0.0, 0.0, D_intra, 13,
                       0.0, EPSILON);
    basis[13] = Sphere(0.5 * a, 0.5 * a, 0.0 * a, radius, 0.0, 0.0, D_intra, 14,
                       0.0, EPSILON);
    /*
            std::vector<Cylinder_XY> basis;
            Lattice lattice(cube_length, cube_length, cube_length, T2_e, 0.0,
       D_extra,100); RPSLatticeInitializer<Cylinder_XY> rpsli(lattice,0);
            rpsli.gammaRadialDist( 12344124, alpha,  beta, .0001,
       cube_length/10); rpsli.uniformCenterDist( 12344213* 5 );
            rpsli.setRegions();
            rpsli.correctEdges();
            lattice = rpsli.lat; //needed to reinitialize basis size (since it
       was initialized to 100 and there will be > 100 cylinders).

            // if (rpsli.basis.size() != 100){std::cout << " Basis Size Does not
       equal 100 " << std::endl;}
            // std::cout << "lattice basis size = " << lattice.getBasisSize() <<
       std::endl; for (int i = 0; i < rpsli.basis.size(); i++){
                    basis.push_back(Cylinder_XY(0.0, 0.0, 0.0,  T2_i,0.0,
       D_intra, i+1, 0.0)); basis[i].setRadius(rpsli.basis[i].getRadius() );
                    basis[i].setCenter(rpsli.basis[i].getCenter() );
                    basis[i].setEPS( (1E-13));
                    basis[i].setRegion(rpsli.basis[i].getRegion());
                    std::cout << rpsli.basis[i].getCenter()  << "  " <<
       rpsli.basis[i].getRadius() << " " << std::endl;

            }

    */

    std::vector<int> plan(3);
    plan[0] = 0;
    plan[1] = NOI;
    plan[2] = NOI;
    std::vector<int> numOfSMPerDevice(1);
    numOfSMPerDevice[0] = 14;
    numOfSMPerDevice[1] = 2;

    pas.runAcquisitionStreamLattice(&basis[0], lattice, timestep, blocks,
                                    threads, 1, plan, numOfSMPerDevice);
    // pas.calcEveryADC();
    // pas.changeSeeds(time(NULL)*(i+1));

    std::cout << std::endl << " Signals " << std::endl;

    for (int i = 0; i < NOI; i++) {
        for (int j = 0; j < NOM; j++) {
            std::cout << setprecision(20);
            std::cout
                << pas.getAcquisition(i).getGradientFunctors()[j].getFreq()
                << " ";
            std::cout << pas.getAcquisition(i).getGradientFunctors()[j].getG()
                      << " ";
            std::cout << pas.getAcquisition(i).getSignal(j) << std::endl;
        }
    }
}
