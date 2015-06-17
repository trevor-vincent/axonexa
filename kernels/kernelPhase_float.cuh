
template <class GradientType, class Basis>
__global__ void updateWalkersPhase(
				 const SimuParamsPhase *par, 
				 const Basis* basis, 
				 const GradientType * G, 
				 curandState* globalState, 
				 real* phase
				 )
{

  const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;

  curandState localState = globalState[tid];
  
  real phi, theta;
  real speed = sqrt(6.0f*basis->getD()/par->timestep);
  real steplength = speed*par->timestep;
	
  #if defined KAHAN_SUMMATION
  real kahanT = 0.0f, kahanC = 0.0f, kahanY = 0.0f, phaseTemp;
  #endif	
  
  Vector3 r = basis[0].unifRand(localState);
  
  for (int i = 0; i < par->measurements; i++){
    //phase[tid + i*par->number_of_particles] = par->phase_initial;
	
	#if defined KAHAN_SUMMATION
	  kahanY = GAMMA*(G[i](0.0f)*r)*par->timestep;
      kahanT = kahanY;
      kahanC = kahanT - kahanY;
      phase[tid + i*par->number_of_particles] = kahanT;
    #else 
	phase[tid + i*par->number_of_particles] = GAMMA*(G[i](0.0f)*r)*par->timestep;
	#endif
  }
	
  for (int i = 1; i < par->steps; i++){
	
    Vector3 ri = r;

#if defined LATTICE_PICKING
	
	r += Vector3(steplength*randomSign(localState),steplength*randomSign(localState),steplength*randomSign(localState));

#else

    phi = 2.0f*PI*curand_uniform( &localState );
    theta = acos(2.0f*curand_uniform( &localState ) - 1);
	
	r += Vector3(speed*par->timestep*sin(theta)*cos(phi), speed*par->timestep*sin(theta)*sin(phi), speed*par->timestep*cos(theta));
	
#endif
   
#if defined SPECULAR_REFLECTION

    boundaryNormal(ri,r, speed, basis, par->timestep);

#else

    if ( !basis[0].inside(r) )
      {
		r = ri;
      }    
    
#endif

    for (int j = 0; j < par->measurements; j++){

	  #if defined KAHAN_SUMMATION
	  phaseTemp = phase[tid + j*par->number_of_particles];
	  kahanY = GAMMA*(G[j](par->timestep*i)*r)*par->timestep - kahanC;
      kahanT = phaseTemp + kahanY;
      kahanC = (kahanT - phaseTemp) - kahanY;
      phase[tid + j*par->number_of_particles] = kahanT;	
	  #else
	  phase[tid + j*par->number_of_particles] += GAMMA*(G[j](par->timestep*i)*r)*par->timestep;
	  #endif
    } 
	
  }

  //globalState[tid] = localState; 

}


