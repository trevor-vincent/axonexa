#define WARP_SIZE 32
#define USE_DOUBLE
//#define SPECULAR_REFLECTION
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
			
		real gradient_duration= 50;	
		real gradient_spacing  = 2.0;
		real echo_time = 2.0*gradient_duration + gradient_spacing ;
		int number_of_timesteps = (int) (echo_time/timestep);		
		phaseAcquisition<CosGFunc> pa(NOM,number_of_timesteps,number_of_particles,j*time(NULL));
		
		for (int i = 0; i < NOM; i++) {
					
			int N = 1+j;
			real G = i*0.0000025*N;
	
			CosGFunc cosGRAD(G, gradient_duration,gradient_spacing, N, Vector3(0.0,0.0,1.0));
			pa.addMeasurement(cosGRAD);	
				
		}
		pas.addAcquisition(pa); 
	}
	
	real radius = .0025;
	real D_extra = 2.5E-6;
	real D_intra = 1.0E-6;
	real T2_i = 200;
	real T2_e = 200;
	// real f = .8;
	// real d = sqrt( 2.0*PI*radius*radius/( sqrt(3.0)*f ) );
	// real a = 0.00198166364880300550667251438256; //sqrt( PI*radius*radius / f );
	
	// Cylinder_XY cyls[4];
	// Lattice lattice(d, 2.0*d*0.86602540378443864676372317075294, d, T2_e, 0.0, D_extra,4);
	// cyls[0] = Cylinder_XY(d/2.0, 0.0,  radius,  T2_i,0.0, D_intra, 1,0.0);
	// cyls[1] = Cylinder_XY(0.0, d*0.86602540378443864676372317075294,  radius,  T2_i,0.0, D_intra, 2,0.0);
	// cyls[2] = Cylinder_XY(d/2.0, 2.0*d*0.86602540378443864676372317075294,  radius,  T2_i,0.0, D_intra, 3,0.0);
	// cyls[3] = Cylinder_XY(d, d*0.86602540378443864676372317075294,  radius,  T2_i,0.0, D_intra, 4,0.0);


	
	// real l = .01;
	// Plane cyls[2];
	// Lattice lattice(l, l, l, T2_e, 0.0, D_extra,2);
	// cyls[0] = Plane(Vector3(0.0,0.0,0.0),Vector3(1.0,0.0,0.0));
	// cyls[1] = Plane(Vector3(l,0.0,0.0),Vector3(-1.0,0.0,0.0));
	
	// Empty nothing[1];
	// Lattice lattice(d, 2.0*d*0.86602540378443864676372317075294, d, T2_e, 0.0, D_extra);
	// nothing[0] = Empty();

	// for (int i = 0; i < NOI; i++){
		// pas.runAcquisitionLattice(i,cyls, lattice, timestep, blocks, threads, 14); 
    // }
	
  std::vector<int> plan(3); plan[0] = 0; plan[1] = NOI;  plan[2] = NOI;
  std::vector<int> numOfSMPerDevice(1); numOfSMPerDevice[0] = 14; numOfSMPerDevice[1] = 2; 
 
 	for (int i = 0; i < 16; i++){
	real f = i*.05;
	std::cout << std::endl << " f = " << f << std::endl << std::endl;
	real a = sqrt( PI*radius*radius / f );
	std::vector<Cylinder_XY> cyls(1);
	Lattice lattice(a, a, a, T2_e, 0.0, D_extra,1);
	cyls[0] = Cylinder_XY(a/2.0, a/2.0,  radius,  T2_i,0.0, D_intra, 1,0.0);
	pas.runAcquisitionLattice(i,&cyls[0], lattice,  timestep, blocks, threads, 14); 
	pas.getAcquisition(i).calcADC();
    std::cout << setprecision(20);
    std::cout << pas.getAcquisition(i).getADC() << std::endl;
	  
	}

 

}
