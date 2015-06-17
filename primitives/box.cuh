//////////////////////////////////////////////////////
// This class only supports rejection sampling
// This class is only used to test against WCDMC
//////////////////////////////////////////////////////

class Box {

private:

	real l_half;
	real w_half;
	real d_half;
	real D;
	
public:

  __device__ __host__ Box(){}
  __device__ __host__ Box(real _l, real _w, real _d, real _D){
	l_half = _l_half;
	w_half = _w_half;
	d_half = _d_half;
	D = _D;
  }
	
  __device__ Vector3 unifRand(curandState & localState) const{
	Vector3 r;
	do{
		r.x = l_half*(2.0f*curand_uniform(&localState) - 1.0f);
		r.y = w_half*(2.0f*curand_uniform(&localState) - 1.0f);
		r.z = d_half*(2.0f*curand_uniform(&localState) - 1.0f);
	} while ( !inside(r) );
	return r;
 }
	
  __device__ bool inside(const Vector3 & r) const{
    return ( fabs(r.x) < l_half &&  fabs(r.y) < h_half && fabs(r.z) < d_half );
  }
  
   __device__ __host__ real getD() const{
		return D;
  }
 
};
