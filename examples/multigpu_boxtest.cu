#define WARP_SIZE 32
#define USE_DOUBLE
// #define SPECULAR_REFLECTION
//#define USE_RELAXATION

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

  int number_of_particles = 57344; //needs to be a factor of two
  real D = 2.5E-6;
  real timestep = .01;  
  real speed = sqrt(6.0*D/timestep);

  //real rad = .005;
  int threads = 128;
  int blocks = number_of_particles/threads;
  
  //for (int z = 2; z < 100; z = z + 2){
  int NOI = 4;
  int NOM = 100;
  phaseAcquisitionStream<CosGFunc> pas(number_of_particles); 
  
	for (int j = 0; j < NOI; j++){
			
		real gradient_duration= 20.0;	
		real gradient_spacing  = 0.0;
		real echo_time = 2.0*gradient_duration + gradient_spacing ;
		int number_of_timesteps = (int) (echo_time/timestep);		
		// std::cout << "number_of_timesteps = " << number_of_timesteps << std::endl;
		phaseAcquisition<CosGFunc> pa(NOM,number_of_timesteps,number_of_particles,j*time(NULL));
		
		for (int i = 0; i < NOM; i++) {
					
			int N = 1+4*j;
			real G = (i+1)*0.0000025*N;
	
			CosGFunc cosGRAD(G, gradient_duration,gradient_spacing, N, Vector3(1.0,0.0,0.0));
			pa.addMeasurement(cosGRAD);	
				
		}
		pas.addAcquisition(pa);	
  }
	
  real rad = .0025;
  Sphere spher(0.0,0.0,0.0,rad,0.0,0.0,D,0,0.0);
  
  std::vector<int> plan(3);
  std::vector<int> numOfSMPerDevice(2); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
 
  // CPUtimer timer1,timer2,timer3,timer4,timer5,timer6;
  // std::cout << " 2 GPU = " << std::endl;  
  // plan[0] = 0; plan[1] = NOI-16;  plan[2] = NOI;
  // timer1.start();
  // pas.runAcquisitionStream(spher, timestep, blocks, threads, 2, plan, numOfSMPerDevice); 
  // timer1.stop();
  // timer1.display();
  // pas.flushADC(); 
  
  plan[0] = 0; plan[1] = NOI;  plan[2] = NOI;
// std::cout << " 1 GPU = " << std::endl;  
  // timer2.start();
  pas.runAcquisitionStream(spher, timestep, blocks, threads, 1, plan, numOfSMPerDevice); 
  // timer2.stop();
  // timer2.display();  
 // pas.flushADC(); 

 
// std::cout << " WC = " << std::endl;
// timer3.start(); 
// for (int i = 0; i < NOI; i++){
	// pas.runAcquisitionWC(i,spher, timestep, blocks, threads);
// } 
// timer3.stop();
// timer3.display();
// pas.flushADC();

// std::cout << " CPU = " << std::endl;
// timer4.start();
// for (int i = 0; i < NOI; i++){
	// pas.runAcquisitionCPU(i,spher, timestep);
// } 
// timer4.stop();
// timer4.display();
// pas.flushADC();
 
  //pas.flushADC();
  // pas.calcEveryADC();
  
  
  
  
  
   // for (int i = 0; i < NOI; i++){
	// std::cout << setprecision(20);
	// std::cout << pas.getAcquisition(i).getGradientFunctors()[0].getFreq() << " " ;
	// std::cout << pas.getAcquisition(i).getADC() << std::endl;
  // }


}
