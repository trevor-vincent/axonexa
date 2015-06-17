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

  int number_of_particles = 14336; //needs to be a factor of two
  real D = 20.0;
  //real D = 1E-6;
  real timestep = .001;
  
  real speed = sqrt(6.0*D/timestep);
  real rad = 50;
  //real rad = .005;
  int threads = 128;
  int blocks = number_of_particles/threads;
	
  magAcquisitionStream<QSC_GRE> pas(number_of_particles);
  // int echo_time = 1;
  // int number_of_timesteps = (int) (echo_time/timestep);
  // real gradient_duration = echo_time*.25;
  // real gradient_spacing = echo_time - 2.0*gradient_duration;	
  // QSC_GRE qsc(.0006, 1.5811388300841896659994467722164e-4,gradient_duration,gradient_spacing);
  // magAcquisition<QSC_GRE> pa(1,number_of_timesteps,number_of_particles,1234);
  // pa.addMeasurement(qsc);
  // pas.addAcquisition(pa); 
  
  int numOfAcq = 20;
  
  for (int i = 0; i < numOfAcq; i++){
	
	  real echo_time = 1;
	  int number_of_timesteps = (int) (echo_time/timestep);
	  real gradient_duration = echo_time*.25;
	  real gradient_spacing = echo_time - 2.0*gradient_duration;	
	  QSC_GRE qsc(0.00000001, i*7.0710678118654752440084436210485e-6/((real)numOfAcq),gradient_duration,gradient_spacing);
	  // std::cout << "echo_time = " << echo_time << std::endl;
	  // std::cout << "gradient_duration = " << gradient_duration << std::endl;
	  // std::cout << "gradient_spacing = " << gradient_spacing << std::endl;
	  
	  //RECT_GRE qsc(0.0, Vector3(i*7.0710678118654752440084436210485e-6/1000.0, 0.0, 0.0),gradient_duration,gradient_spacing);
	  magAcquisition<QSC_GRE> pa(1,number_of_timesteps,i*time(NULL));
	  pa.addMeasurement(qsc);
	  pas.addAcquisition(pa); 
	
  } 
  
  //pas.runCPUAcquisition(0, rad, timestep,  c, speed);	
  Cylinder_XY cylinder(0.0,0.0,rad,0.0,0.0,D,0);
  Vector3 initialM(1.0,0.0,0.0);
  Substrate sub(90,90,.001);
  // std::vector<int> plan(2); plan[0] = 0; plan[1] = 4; 
  // std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; 
  
  std::vector<int> plan(3); plan[0] = 0; plan[1] = 20;  plan[2] = 20;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
  
  pas.runAcquisitionStream(cylinder,  sub, initialM,  timestep, blocks, threads, 1, plan, numOfSMPerDevice);

	CPUtimer timer;
	timer.start();
	for (int i = 0; i < numOfAcq; i++){
	  pas.runCPUAcquisition(i, cylinder,  sub, initialM,  timestep);
  }
  timer.stop();
  timer.display();
  
  //pas.runAcquisition(0,cylinder,  sub, initialM,  timestep, blocks, threads,1,2);
  // pas.runAcquisition(1,cylinder,  sub, initialM,  timestep, blocks, threads);
   // pas.runAcquisition(2,cylinder,  sub, initialM,  timestep, blocks, threads);
    // pas.runAcquisition(3,cylinder,  sub, initialM,  timestep, blocks, threads);
  
  //pas.runCPUAcquisition(0,cylinder,  sub, initialM,  timestep);
  //pas.runCPUAcquisition(0,cylinder,  sub, initialM,  timestep);
  // pas.getAcquisition(0).calcTotalMagnetization(number_of_particles);
  // cout << pas.getAcquisition(0).getTotalMagnetization()[0].magnitude() << endl;
  //pas.runCPUAcquisition(0, rad, timestep,  c, speed);
  //pas.runStreamOfAcquisitions(rad, timestep, c, speed,  blocks,  threads);
  //pas.calcEveryTotalMag();
  
  //pas.calcEveryTotalMag();
  // for (int i = 0; i < numOfAcq; i++){
  // cout << setprecision(20) << endl;
  // cout << pas.getAcquisition(i).getFieldFunctors()[0].G<< " " ;
  // cout << Vector3( pas.getAcquisition(i).getMx()[0], pas.getAcquisition(i).getMy()[0], pas.getAcquisition(i).getMz()[0])*(1.0/number_of_particles);
  // }
  
  
  
  
  
  
  // for (int i = 0; i < 200; i++){
  // cout << i*1.5811388300841896659994467722164e-7*i*1.5811388300841896659994467722164e-7 << " " << pas.acqSet[i].totalMagnetization[0].magnitude() << endl;
  // }
	
  //pas.runAcquisition(0, number_of_particles, rad, timestep, c, speed,  blocks,  threads);
  //pas.runAcquisition(1, number_of_particles, rad, timestep, c, speed,  blocks,  threads);
  //pas.runAcquisition(0, number_of_particles, rad, timestep, c, speed,  blocks,  threads, 1231222393);
	
  //ms.addMeasurementSet(number_of_timesteps,2);
  //ms.addMeasurement(rg1,timestep); ms.addMeasurement(rg2,timestep);
  //ms.runStreamOfAcquisitions(number_of_particles, rad, timestep,  c, speed,  blocks,  threads, 12343432);
  //ms.runAcquisition(0, number_of_particles, rad, timestep, c, speed,  blocks,  threads);
  //ms.runAcquisition(1, number_of_particles, rad, timestep, c, speed,  blocks,  threads);
}
