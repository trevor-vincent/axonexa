//get rid of substrate

#include <assert.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <thrust/sort.h>
#include <cmath>
#include <thrust/device_ptr.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>

#define USE_DOUBLE
//#define KERNEL_DEBUG

#define WARP_SIZE 32
#if defined USE_DOUBLE
typedef double real;
#define EPSILON 1e-14

#else
typedef float real;
#define EPSILON 1e-6

#endif

using namespace std;


//#define SPECULAR_REFLECTION
#define GAMMA 267500.0 // ms^-1 * T^-1
#define PI 3.1415926535897932384626433832795
//#define SPECULAR_REFLECTION


#include "misc.cuh"
#include "vector3.cuh"
#include "cudaVector.cuh"
#include "timer.cuh"
#include "compare.cuh"
#include "pinnedVector.cuh"
#include "cudaVector.cu"
#include "pinnedVector.cu"
#include "bfunctors.cuh"

//#include "intersect.cuh"

#include "substrate.cuh"
#include "cylinderXY.cuh"
#include "Sphere.cuh"
#include "simuparams.cuh"
#include "kernelMag.cuh"
#include "kernelDEBUG.cuh"
#include "kernelPhase.cuh"
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
  real D = 2.5E-6;
  //real D = 1E-6;
  real timestep = .001;  
  real speed = sqrt(6.0*D/timestep);

  //real rad = .005;
  int threads = 128;
  int blocks = number_of_particles/threads;
	
  std::cout << "NOI" << " " << " TDMC " << " " << " WCDMC " << std::endl;
	
  for (int z = 0; z < 100; z += 5){
  
  int NOI = z+1;
  int NOM = 2;
  phaseAcquisitionStream<CosGFunc> pas(number_of_particles); 
  
	for (int j = 0; j < NOI; j++){
			
		real gradient_duration= 10;	
		real gradient_spacing  = 2.0;
		real echo_time = 2.0*gradient_duration + gradient_spacing ;
		int number_of_timesteps = (int) (echo_time/timestep);		
		phaseAcquisition<CosGFunc> pa(NOM,number_of_timesteps,number_of_particles,j*time(NULL));
		
		for (int i = 0; i < NOM; i++) {
					
			int N = 1+4*j;
			real G = i*0.0000025*N;
	
			CosGFunc cosGRAD(G, gradient_duration,gradient_spacing, N, Vector3(1.0,0.0,0.0));
			pa.addMeasurement(cosGRAD);	
				
		}
		pas.addAcquisition(pa); 
		
	}

  Sphere spher(0.0,0.0,0.0,rad,0.0,0.0,D,0);
  
  std::vector<int> plan(3); plan[0] = 0; plan[1] = NOI;  plan[2] = NOI;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
 
CPUtimer timer1, timer2;
timer1.start();
pas.runAcquisitionStream(spher, timestep, blocks, threads, 1, plan, numOfSMPerDevice); 
timer1.stop();

timer2.start(); 
for (int i = 0; i < NOI; i++){
	pas.runAcquisitionWC(i,spher, timestep, blocks, threads);
} 
timer2.stop();

std::cout << NOI << " " << timer1.getTime() << " " << timer2.getTime() << std::endl;
  
}

  
  
}
