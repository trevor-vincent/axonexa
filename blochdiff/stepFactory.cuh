class fixedSphereStep {
	__device__ __host__ void update(Vector3 & r, curandState & localState){
		phi = 2.0*PI*curand_uniform( &localState );
		theta = acos(2.0*curand_uniform( &localState ) - 1);
		r += Vector3(speed*par->timestep*sin(theta)*cos(phi), 
			         speed*par->timestep*sin(theta)*sin(phi), 
					 speed*par->timestep*cos(theta));
	}
};

class specularSinglePore {
	__device__ __host__ void amend(Vector3 & r, curandState & localState){
		phi = 2.0*PI*curand_uniform( &localState );
		theta = acos(2.0*curand_uniform( &localState ) - 1);
		r += Vector3(speed*par->timestep*sin(theta)*cos(phi), 
			         speed*par->timestep*sin(theta)*sin(phi), 
					 speed*par->timestep*cos(theta));
	}




}