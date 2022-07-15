#define WARP_SIZE 32
#define USE_DOUBLE
#define SPECULAR_REFLECTION
//#define USE_RELAXATION
#define GAMMA 267500.0 // ms^-1 * T^-1
//#define GAMMA (-203789.0)
#define PI 3.1415926535897932384626433832795
#define D_NOT_ZERO
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
//#define EPSILON 1e-14 //good for small radii
#define EPSILON 1e-14 //good for large radii

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
#include "calcf.cuh"


int main (){
 
  cudaFuncSetCacheConfig( "updateWalkersMag", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "setup_kernel", cudaFuncCachePreferL1 );
  cudaFuncSetCacheConfig( "_functionReduceAtom", cudaFuncCachePreferShared );

  int number_of_particles = 114688; //needs to be a factor of two
  real D = 2.5E-6;
  //real D = 1E-6;
  real timestep = .0001;
  
  real speed = sqrt(6.0*D/timestep);
  real rad = .0025;
  //real rad = .005;
  int threads = 128;
  int blocks = number_of_particles/threads;
	
  // int echo_time = 1;
  // int number_of_timesteps = (int) (echo_time/timestep);
  // real gradient_duration = echo_time*.25;
  // real gradient_spacing = echo_time - 2.0*gradient_duration;	
  // QSC_GRE qsc(.0006, 1.5811388300841896659994467722164e-4,gradient_duration,gradient_spacing);
  // magAcquisition<QSC_GRE> pa(1,number_of_timesteps,number_of_particles,1234);
  // pa.addMeasurement(qsc);
  // pas.addAcquisition(pa); 
  
  int numOfAcq = 1000;
  std::vector<real> B0v(1);
   // B0v[0] = .00001;
   // B0v[1] = .00005;
   // B0v[2] = .0001;
   // B0v[3] = .0005;
   // B0v[4] = .001;
   // B0v[5] = .005;
   // B0v[6] = .01;
   // B0v[7] = .05;
   // B0v[8] = .1;
   // B0v[9] = .5;
   
   //B0v[0] = .00005;
   //B0v[1] = .0001;
   B0v[0] = .0005;
 
  // for (int r = 0; r < 10; r++){
  for (int r = 0; r < 1; r++){
  real B0 = B0v[r]; 
  std::cout << std::endl << " B0 = " << B0 << std::endl << std::endl;
  int repeats = 3;

  magAcquisitionStream<QSC_GRE> pas(number_of_particles);
  real maxDuration = 0.05;

  for (int i = 0; i < numOfAcq; i++){

	  // real gradient_duration = (i+1)*maxDuration/numOfAcq; 
	  real gradient_duration = (i+1)*maxDuration/numOfAcq; 
      real gradient_spacing = 0.01; 
	  real echo_time = 2.0*gradient_duration + gradient_spacing;
	  real G1 = 0.0;
	  real G2 = .000006;
	  int number_of_timesteps = (int) (echo_time/timestep);
	  QSC_GRE qsc1(B0,G1,gradient_duration,gradient_spacing);
	  QSC_GRE qsc2(B0,G2,gradient_duration,gradient_spacing);
	  
	  magAcquisition<QSC_GRE> pa(2,number_of_timesteps,i*time(NULL));
	  pa.addMeasurement(qsc1);
	  pa.addMeasurement(qsc2);
	  
	  pas.addAcquisition(pa); 
	
  } 
  
  //pas.runCPUAcquisition(0, rad, timestep,  c, speed);	
  Sphere spher(0.0,0.0,0.0,rad,0.0,0.0,D,0,0.0);
  Vector3 initialM(1.0,0.0,0.0);
  
  std::vector<int> plan(3); plan[0] = 0; plan[1] = numOfAcq;  plan[2] = numOfAcq;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
  
  for (int i = 0; i < repeats; i++){
	pas.runAcquisitionStream(spher, initialM,  timestep, blocks, threads, 1, plan, numOfSMPerDevice);
	pas.calcEveryADC();
	pas.changeSeeds(time(NULL)*(i+1));
  }
  
 // for (int i = 0; i < numOfAcq; i++){
	// pas.runCPUAcquisition(i,cylinder, initialM,  timestep);
 // }
 
 for (int i = 0; i < numOfAcq; i++){
	std::cout << setprecision(20);
	std::cout << pas.getAcquisition(i).getFieldFunctors()[0].getDiffusionTime() << " " ;
	std::cout << pas.getAcquisition(i).getADC() << " " ;
	std::cout << pas.getAcquisition(i).getADCError() << std::endl;
  }
  
  }
}
