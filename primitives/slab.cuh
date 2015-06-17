//////////////////////////////////////////////////////
// This class only supports rejection sampling
// This class is only used to test against WCDMC
//////////////////////////////////////////////////////

class Slab {

private:

	real l; //length
	real D; //Diffusion coefficient
	
public:

  __device__ __host__ Slab(){}
  __device__ __host__ Slab(real _l, real _D){
	l = _l;
	D = _D;
  }

  __device__ Vector3 unifRand(curandState & localState) const{
	Vector3 r;
	do{
		r.x = l*(2.0f*curand_uniform(&localState) - 1.0f);
		r.y = l*(2.0f*curand_uniform(&localState) - 1.0f);
		r.z = l*(2.0f*curand_uniform(&localState) - 1.0f);
	} while ( !inside(r) );
	return r;
   }
 
  __host__ Vector3 unifRandCPU() const{
	Vector3 r;
	do{
		r.x = l*(2.0f*unifRandCPP()  - 1.0f);
		r.y = l*(2.0f*unifRandCPP()  - 1.0f);
		r.z = l*(2.0f*unifRandCPP()  - 1.0f);
	} while ( !inside(r) );
	return r;
  }	
  
  __device__ __host__ bool inside(const Vector3 & r) const{
    return ( r.x < l && r.x > 0 );
  }
  
 __device__ __host__ bool inside(const real x, const real y, const real z) const{
    return ( x < l && x > 0 );
  }
  
   __device__ __host__ real getD() const{
		return D;
  }
 
};
