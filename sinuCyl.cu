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
#define EPSILON 5e-12

#else
typedef float real;
#define EPSILON 1e-6

#endif

using namespace std;


#define SPECULAR_REFLECTION
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
#include "simuparams.cuh"
#include "kernels.cuh"
#include "kernelDEBUG.cuh"
#include "CPUkernels.cuh"
#include "magAcquisition.cuh"
#include "magAcquisitionStream.cuh"
#include "blochdiff.cuh"



int main (){

	cudaFuncSetCacheConfig( "updateWalkersMag", cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( "setup_kernel", cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( "_functionReduceAtom", cudaFuncCachePreferShared );

  int number_of_particles = 14336; //needs to be a factor of two
  real D = 2.5E-6;
  real timestep = .001;
  real speed = sqrt(6.0*D/timestep);
  real rad = .001;
  int threads = 128;
  int blocks = number_of_particles/threads;
	
  magAcquisitionStream<COS_GRE> pas(number_of_particles);
  
  int NOI = 4;
  int NOM = 2;

	for (int j = 0; j < NOI; j++){
			
		real gradient_duration= 4;	
		real gradient_spacing  = 2.0;
		real echo_time = 2.0*gradient_duration + gradient_spacing ;
		int number_of_timesteps = (int) (echo_time/timestep);		
		magAcquisition<COS_GRE> pa(NOM,number_of_timesteps,j*time(NULL));
		
		for (int i = 0; i < NOM; i++) {
					
			int N = 1+4*j;
			real G = i*0.0000025*N;
	
			COS_GRE cosGRAD(0.1, G, Vector3(1.0,0.0,0.0), gradient_duration,gradient_spacing, N);
			pa.addMeasurement(cosGRAD);	
				
		}
		pas.addAcquisition(pa); 
	}

	
  Cylinder_XY cylinder(0.0,0.0,rad,0.0,0.0,D,0);
  Vector3 initialM(1.0,0.0,0.0);
  Substrate sub(.001,.001,.001);
  
  std::vector<int> plan(3); plan[0] = 0; plan[1] = 2;  plan[2] = 4;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
  
 pas.runAcquisitionStream(cylinder,  sub, initialM,  timestep, blocks, threads, 2, plan, numOfSMPerDevice); 

  // pas.runAcquisition(0, cylinder,  sub, initialM,  timestep, blocks, threads, 1 , 2);
  // pas.runAcquisition(1, cylinder,  sub, initialM,  timestep, blocks, threads, 1 , 2);
  // pas.runAcquisition(2, cylinder,  sub, initialM,  timestep, blocks, threads, 0 , 14);
  // pas.runAcquisition(3, cylinder,  sub, initialM,  timestep, blocks, threads, 0 , 14);
  
  for (int i = 0; i < NOI; i++){
  // std::cout << " ERROR HERE " << std::endl;
    // cout << Vector3( pas.getAcquisition(i).getMx()[0], pas.getAcquisition(i).getMy()[0], pas.getAcquisition(i).getMz()[0])*(1.0/number_of_particles);
	pas.getAcquisition(i).calcADC();
	std::cout << setprecision(20);
	std::cout << pas.getAcquisition(i).getFieldFunctors()[0].getFreq() << " " ;
	std::cout << pas.getAcquisition(i).getADC() << std::endl;
  }
  
  
}
