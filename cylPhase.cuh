#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
//#define USE_RELAXATION

#include <assert.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>

#if defined USE_DOUBLE
typedef double real;
#define EPSILON 1e-14

#else
typedef float real;
#define EPSILON 1e-6

#endif

using namespace std;

#include "misc.cuh"
#include "vector3.cuh"
#include "cudaVector.cuh"
#include "timer.cuh"
#include "compare.cuh"
#include "pinnedVector.cuh"
#include "cudaVector.cu"
#include "pinnedVector.cu"
#include "bfunctors.cuh"
#include "substrate.cuh"
#include "cylinderXY.cuh"
#include "Sphere.cuh"
#include "plane.cuh"
#include "empty.cuh"
#include "lattice.cuh"
#include "simuparams.cuh"
#include "boundaryCheck.cuh"
#include "kernelSetup.cuh"
#include "kernelMag.cuh"
#include "kernelDEBUG.cuh"
#include "kernelPhase.cuh"
#include "kernelLattice.cuh"
#include "kernelWC.cuh"
#include "CPUkernels.cuh"
#include "gfunctors.cuh"
#include "phaseAcquisition.cuh"
#include "phaseAcquisitionStream.cuh"
#include "magAcquisition.cuh"
#include "magAcquisitionStream.cuh"
#include "blochdiff.cuh"



int main (){

	cudaFuncSetCacheConfig( "updateWalkersMag", cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( "setup_kernel", cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( "_functionReduceAtom", cudaFuncCachePreferShared );

  int number_of_particles = 114688; //needs to be a factor of two
  real timestep = .001;  

  int threads = 128;
  int blocks = number_of_particles/threads;
	
  phaseAcquisitionStream<CosGFunc> pas(number_of_particles);
  
   int NOI = 1;
   int NOM = 2;

	for (int j = 0; j < NOI; j++){
			
		real gradient_duration= 10;	
		real gradient_spacing  = 2.0;
		real echo_time = 2.0*gradient_duration + gradient_spacing ;
		int number_of_timesteps = (int) (echo_time/timestep);		
		phaseAcquisition<CosGFunc> pa(NOM,number_of_timesteps,number_of_particles,j*time(NULL));
		
		for (int i = 0; i < NOM; i++) {
					
			int N = 1+j;
			real G = i*0.0000025*N;
	
			CosGFunc cosGRAD(G, gradient_duration,gradient_spacing, N, Vector3(1.0,0.0,0.0));
			pa.addMeasurement(cosGRAD);	
				
		}
		pas.addAcquisition(pa); 
	}
	
	std::vector<real> lengths(4);
	lengths[0] = .001;
	lengths[1] = .002;
	lengths[2] = .005;
	lengths[3] = .010;
	for (int r = o; r < 4; r++){
	
	real D_extra = 2.5E-6;
	real D_intra = 1.0E-6;
	real T2_i = 200;
	real T2_e = 200;
	
	real l = lengths[r];
	std::cout << " CYLINDER DIAMETER = " << l << std::endl;
    Cylinder_XY cyls(0.0,0.0,l/2.0,0.0,0.0,D,0,0.0);
	
	// Empty nothing[1];
	// Lattice lattice(d, 2.0*d*0.86602540378443864676372317075294, d, T2_e, 0.0, D_extra);
	// nothing[0] = Empty();

	// for (int i = 0; i < NOI; i++){
		// pas.runAcquisitionLattice(i,cyls, lattice, timestep, blocks, threads, 14); 
    // }
	
  std::vector<int> plan(3); plan[0] = 0; plan[1] = NOI;  plan[2] = NOI;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
 
  int repeats = 3;
 
  for (int i = 0; i < repeats; i++){
	  pas.runAcquisitionStream(cyls, timestep, blocks, threads, 1, plan, numOfSMPerDevice); 
	  pas.calcEveryADC();
	  pas.changeSeeds(time(NULL)*(i+1));
  }

 
 for (int i = 0; i < NOI; i++){
	std::cout << setprecision(20);
	std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq() << " " ;
	std::cout << pas.getAcquisition(i).getADC() << std::endl;
  }
  
 }
}
