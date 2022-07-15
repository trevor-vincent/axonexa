#define WARP_SIZE 32
#define KAHAN_SUMMATION
//#define USE_DOUBLE
//#define SPECULAR_REFLECTION
//#define USE_RELAXATION
//#define GAMMA 267500.0 // ms^-1 * T^-1
//#define PI 3.1415926535897932384626433832795
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
#define GAMMA 267500.0
#define PI 3.1415926535897932384626433832795
#else
typedef float real;
#define EPSILON 1e-6
#define GAMMA 267500.0f
#define PI 3.1415926535897932384626433832795f
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
#include "slab.cuh"
#include "lattice.cuh"
#include "simuparams.cuh"

#if defined USE_DOUBLE
#include "boundaryCheck.cuh"
#else 
#include "boundaryCheck_float.cuh"
#endif

#include "kernelSetup.cuh"
#include "kernelMag.cuh"
#include "kernelDEBUG.cuh"

#if defined USE_DOUBLE
#include "kernelPhase.cuh"
#else
#include "kernelPhase_float.cuh"
#endif

#include "kernelLattice.cuh"
#include "kernelWC.cuh"
#include "CPUkernels.cuh"
#if defined USE_DOUBLE
#include "gfunctors.cuh"
#else
#include "gfunctors_float.cuh"
#endif

#include "phaseAcquisition.cuh"
#include "phaseAcquisitionStream.cuh"
#include "magAcquisition.cuh"
#include "magAcquisitionStream.cuh"

#if defined USE_DOUBLE
#include "blochdiff.cuh"
#else
#include "blochdiff_float.cuh"
#endif

int main (){

  cudaFuncSetCacheConfig( "updateWalkersLattice", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "updateWalkersPhase", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "updateWalkersMag", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "setup_kernel", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "_functionReduceAtom", cudaFuncCachePreferShared );
  cudaFuncSetCacheConfig( "_functionTransformAndReduceAtom", cudaFuncCachePreferShared );
  cudaFuncSetCacheConfig( "_functionTransformAndSumTwoVectorsAtom", cudaFuncCachePreferShared );
	
  int number_of_particles = 65536; //needs to be a factor of two
  real D = 1.0E-6;
  real timestep = .001;  

  //real rad = .005;
  int threads = 128;
  int blocks = number_of_particles/threads;
  

  phaseAcquisitionStream<SinGFunc> pas(number_of_particles);
  
   int NOI = 100;
   int NOM = 1;

 	real gradient_duration = 10;	
	real gradient_spacing  = 1.0;
	real echo_time = 2.0*gradient_duration + gradient_spacing ;
	int number_of_timesteps = (int) (echo_time/timestep) + 1;		
    phaseAcquisition<SinGFunc> pa(NOM*NOI,number_of_timesteps,number_of_particles,time(NULL));
	real G = .01;
	
	for (int j = 0; j < NOI; j++){
		for (int i = 0; i < NOM; i++) {
			int N = j + 1;
			SinGFunc sinGRAD(G, gradient_duration, gradient_spacing, N, Vector3(1.0,0.0,0.0));
			pa.addMeasurement(sinGRAD);			
		}
	}
	
  pas.addAcquisition(pa);
	
 
  Slab slabby(.005,D);

  for (int i = 0; i < 3; i++){  

  CPUtimer timer1, timer2, timer3;
  timer1.start();
  pas.runAcquisition(0, slabby, timestep, blocks, threads, 14); 
  timer1.stop();
  timer1.display();

  timer2.start(); 
  pas.runAcquisitionWC(0, slabby, timestep, blocks, threads);
  timer2.stop(); 
  timer2.display();
  
  // timer3.start(); 
  // pas.runAcquisitionCPU(0, slabby, timestep);
  // timer3.stop(); 
  // timer3.display();
  
  }
  
  
}
