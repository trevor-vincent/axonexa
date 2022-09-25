#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
//#define USE_INITIALIZE_IN_REGION 0
//#define USE_RELAXATION

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cuda.h>
#include <curand_kernel.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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

void getNGData(std::vector<real> &N, std::vector<real> &G, std::string file) {

    std::ifstream ifs(file.c_str());
    std::string temp;
    real Gt, Nt;
    real gmax = 4.30357 * 1000. / 4.26;
    while (getline(ifs, temp)) {
        std::stringstream ss(temp);
        ss >> Nt >> Gt;
        Gt /= 100.;
        Gt *= gmax;
        Gt *= 1E-6;
        N.push_back(Nt);
        G.push_back(Gt);
    }
}

int main() {

    cudaFuncSetCacheConfig("updateWalkersMag", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("setup_kernel", cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig("_functionReduceAtom", cudaFuncCachePreferShared);

    int number_of_particles = 57344; // needs to be a factor of two
    real timestep = .001;

    int threads = 128;
    int blocks = number_of_particles / threads;

    phaseAcquisitionStream<SinGFunc> pas(number_of_particles);

    std::vector<real> numOfPeriods;
    std::vector<real> gradStrengths;

    getNGData(numOfPeriods, gradStrengths, "data.txt");

    real gradient_duration = 20;
    real gradient_spacing = 4.52;
    real echo_time = 2.0 * gradient_duration + gradient_spacing;
    int number_of_timesteps = (int)(echo_time / timestep);

    phaseAcquisition<SinGFunc> pa(numOfPeriods.size(), number_of_timesteps,
                                  number_of_particles, time(NULL));

    for (int j = 0; j < numOfPeriods.size(); j++) {
        int N = numOfPeriods[j];
        real G = gradStrengths[j];
        std::cout << N << " " << G << std::endl;
        SinGFunc cosGRAD(G, gradient_duration, gradient_spacing, N,
                         Vector3(1.0, 0.0, 0.0));
        pa.addMeasurement(cosGRAD);
    }

    pas.addAcquisition(pa);

    phaseAcquisitionStream<SinGFunc> pas1, pas2;
    pas1 = pas;
    pas1.getAcquisition(0).getSeed() *= 2;
    pas2 = pas;
    pas2.getAcquisition(1).getSeed() *= 3;

    real radius = .0015;
    real D_extra = 2.5E-6;
    real D_intra = 1.0E-6;
    real f = .4;
    real a = pow((16. / 3.) * PI * radius * radius * radius / f, 1. / 3.);

    std::vector<Sphere> basis(1);
    Lattice lattice(a, a, a, 0.0, 0.0, D_extra, 1);
    basis[0] = Sphere(40000.0 * a, 40000.0 * a, 40000.0 * a, 0.0, 0.0, 0.0, 0.0,
                      1, 0.0, EPSILON);

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

    // std::vector<int> plan(3); plan[0] = 0; plan[1] = 1;  plan[2] = 1;
    // std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14;
    // numOfSMPerDevice[1] = 2;

    // pas.runAcquisitionStreamLattice(&basis[0], lattice,  timestep, blocks,
    // threads, 1, plan, numOfSMPerDevice);
    pas.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                              14);
    pas1.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);
    pas2.runAcquisitionLattice(0, &basis[0], lattice, timestep, blocks, threads,
                               14);

    std::cout << std::endl << " Signals " << std::endl;

    for (int j = 0; j < numOfPeriods.size(); j++) {
        std::cout << setprecision(20);
        std::cout << pas.getAcquisition(0).getGradientFunctors()[j].getFreq()
                  << " ";
        std::cout << pas.getAcquisition(0).getGradientFunctors()[j].getG()
                  << " ";
        std::cout << pas.getAcquisition(0).getMx()[j] << " "
                  << pas.getAcquisition(0).getMy()[j] << " "
                  << pas1.getAcquisition(0).getMx()[j] << " "
                  << pas1.getAcquisition(0).getMy()[j] << " "
                  << pas2.getAcquisition(0).getMx()[j] << " "
                  << pas2.getAcquisition(0).getMy()[j] << std::endl;
    }
}
