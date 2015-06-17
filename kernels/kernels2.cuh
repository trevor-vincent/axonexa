//returns the final magnetization vector (so use = in the main code with this result)
/*
  The long drawn out way
  // real phi = GAMMA*(B.magnitude())*timestep;
  // Vector3 C = B*(M*B);
  // Vector3 u = M - C;
  // Vector3 v = B % u;
  // Vector3 r = u*cos(-phi) + v*sin(-phi);
  // return (C + r);

  */

__device__ __host__ Vector3 bloch_noT2(Vector3 M, Vector3 B){

  return (M%B)*GAMMA;

}


__device__ __host__ Vector3 updateMag_noT2_rotate(Vector3 M, Vector3 B, real timestep){
 
  real Bmag = B.magnitude();
  Vector3 C = (B*(1.0/Bmag))*(M*(B*(1.0/Bmag)));
  Vector3 u = M-C;
  return (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
}

__device__ __host__ Vector3 updateMag_noT2_rotate(real* Mx, real* My, real* Mz, Vector3 B, real timestep){
  Vector3 M(*Mx, *My, *Mz);
  real Bmag = B.magnitude();
  Vector3 C = (B*(1.0/Bmag))*(M*(B*(1.0/Bmag)));
  Vector3 u = M-C;
  Vector3 finalM =  (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
  *Mx = finalM.x;
  *My = finalM.y;
  *Mz = finalM.z;
}

template <class FieldType>
__device__ __host__ Vector3 updateMag_noT2_RK4(real* Mx, 
					       real* My, 
					       real* Mz, 
					       FieldType B, 
					       real i,
					       real h)
{

  real tn = i*timestep;
  Vector3 M(*Mx,*My,*Mz);
  Vector3 k1 = (M%B)*timestep*GAMMA;
  Vector3 k2 = bloch_noT2((M + k1*.5), B(tn + .5*h));
  Vector3 k3 = bloch_noT2((M + k2*.5), B(tn + .5*h));
  Vector3 k4 = bloch_noT2((M + k3),B(tn+h));
  Vector3 finalM = M + (k1 + k2*2.0 + k3*2.0 + k4)*(1.0/6.0);

}
  

//relies on the fact that the Magnetization is initialized as Vector3(1.0,0.0,0.0)
//also relies on the fact that there is a Bfield functor in cosntant memory
template <class FieldType, class Basis>
__global__ void updateWalkersMag(const int number_of_particles,
				 const real timestep,  
				 const int number_of_timesteps, 
				 const real c, 
				 const real speed, 
				 const Basis* basis, 
				 Vector3* M, 
				 const FieldType * B, 
				 curandState* globalState, 
				 const int num_of_meas )
{

  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
  curandState localState = globalState[tid];
  //real dt;
  real phi, theta;
	
  Vector3 x = basis.unifRand(&localstate, c);
  
  for (int i = 0; i < num_of_meas; i++){
    M[tid + i*number_of_particles] = Vector3(1.0,0.0,0.0);
    M[tid + i*number_of_particles] = updateMag_noT2_rotate( M[tid + i*number_of_particles], B[i](x,0.0), timestep );
  }
	
  for (int i = 1; i < number_of_timesteps; i++){
	
    Vector3 xi = x; //dt = 0.0;
	
    phi = 2.0*PI*curand_uniform( &localState );
    theta = acos(2.0*curand_uniform( &localState ) - 1);
	
    x += Vector3(speed*timestep*sin(theta)*cos(phi), speed*timestep*sin(theta)*sin(phi), speed*timestep*cos(theta));
    if (!basis.inside(x)){x = xi;}

    for (int j = 0; j < num_of_meas; j++){
      M[tid + j*number_of_particles] = updateMag_noT2_rotate( M[tid + j*number_of_particles], B[j](x,i*timestep) , timestep);
    } 
	
  }
  globalState[tid] = localState; 
}


template <class FieldType, class Basis>
__global__ void updateWalkersMag(const int number_of_particles, 
				 const real timestep,  
				 const int number_of_timesteps, 
				 const real c, 
				 const real speed, 
				 const Basis* basis, 
				 real* Mx, 
				 real* My, 
				 real* Mz,
				 const FieldType * B, 
				 curandState* globalState, 
				 const int num_of_meas ){

  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
  curandState localState = globalState[tid];
  //real dt;
  real phi, theta;
	
  Vector3 x = basis.unifRand(&localstate, c);
  
  for (int i = 0; i < num_of_meas; i++){
    updateMag_noT2_rotate(  &Mx[tid + i*number_of_particles] ,  &My[tid + i*number_of_particles] , &Mz[tid + i*number_of_particles] , B[i](x,0.0), timestep );
  }
	
  for (int i = 1; i < number_of_timesteps; i++){
	
    Vector3 xi = x; //dt = 0.0;
	
    phi = 2.0*PI*curand_uniform( &localState );
    theta = acos(2.0*curand_uniform( &localState ) - 1);
	
    x += Vector3(speed*timestep*sin(theta)*cos(phi), speed*timestep*sin(theta)*sin(phi), speed*timestep*cos(theta));
    if (!basis.inside(x)){x = xi;}

    for (int j = 0; j < num_of_meas; j++){
      updateMag_noT2_rotate(  &Mx[tid + j*number_of_particles] ,  &My[tid + j*number_of_particles] , &Mz[tid + j*number_of_particles] ,B[j](x,i*timestep), timestep );

    } 
	
  }
  globalState[tid] = localState; 
}




/*

  Single timestep kernel

*/


template <class Basis>
__global__ void initWalkers(real *x,
			    real *y,
			    real *z,
			    curandState* globaState,
			    const Basis* basis)
{

  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
   curandState localState = globalState[tid];
   Vector3 r = basis.unifRand(&localstate, c);
   x[tid] = r.x;
   y[tid] = r.y;
   z[tid] = r.z;
}
			    
__global__ void initMagnetization(real *mx,
			    real *my,
			    real *mz,
			    const Vector3 m)
{

  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
   mx[tid] = m.x;
   my[tid] = m.y;
   mz[tid] = m.z;
}


template <class FieldType, class Basis>
__global__ void updateWalkersMag(const int number_of_particles, 
				 const real timestep,  
				 const real c, 
				 const real speed, 
				 const Basis* basis, 
				 real *x,
				 real *y,
				 real *z,
				 real* Mx, 
				 real* My, 
				 real* Mz, 
				 const FieldType * B, 
				 curandState* globalState, 
				 const int num_of_meas ){

  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
  curandState localState = globalState[tid];

  Vector3 x(x[tid], y[tid], z[tid]);
  Vector3 xi = x; //dt = 0.0;
	
  phi = 2.0*PI*curand_uniform( &localState );
  theta = acos(2.0*curand_uniform( &localState ) - 1);
	
  x += Vector3(speed*timestep*sin(theta)*cos(phi), speed*timestep*sin(theta)*sin(phi), speed*timestep*cos(theta));
  if (!basis.inside(x)){x = xi;}

  for (int j = 0; j < num_of_meas; j++){
    updateMag_noT2_rotate(  &Mx[tid + j*number_of_particles] ,  &My[tid + j*number_of_particles] , &Mz[tid + j*number_of_particles] ,B[j](x,i*timestep), timestep );
  } 
	
  globalState[tid] = localState; 
}
