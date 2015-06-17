template < 
		 class FieldType, 
		 class Basis, 
		 class Integrator,
		 class StepCreator,
		 class StepAmender,
	     >
		
__global__ void updateWalkersMag(
					const Initializer *initializer, 
					const Basis* basis, 
					const FieldType * B, 
					curandState* globalState, 
					real* Mx, 
					real* My, 
					real* Mz	 
				 )
{

  Integrator integrator; 
  StepCreator stepCreator; 
  StepAmender stepAmender;
  
  const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
  curandState localState = globalState[tid];
  
  Vector3 r;
  real v;
  
  initializer.initialize(r,v,basis);
	
  for (int i = 1; i < par->steps; i++){
	
    Vector3 ri = r;
	stepCreator.update(r,v,basis);
	stepAmender.correct(ri,r,v,basis,par->timestep);

    for (int j = 0; j < par->measurements; j++){
	integrator.integrate (  
				  &Mx[tid + j*par->number_of_particles],  
			      &My[tid + j*par->number_of_particles], 
			      &Mz[tid + j*par->number_of_particles], 
				  basis,
				  B[j],
			      r,
				  i,
			      par->timestep 
				)	
    } 
	
  }

}