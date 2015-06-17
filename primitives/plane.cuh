class Plane {

private:

  Vector3 center;
  Vector3 normal;
  real permeability;

public:

  __device__ __host__ Plane(){}
  __device__ __host__ Plane(Vector3 _center, Vector3 _normal){
		
		center = _center;
		normal = _normal;
	
	}
  
  __device__ __host__ ~Plane(){}

		
  __device__ __host__ bool intersect(const Vector3 & ri, const Vector3 & rf, real & v) const
  {
		//intersection with entire plane defined by triangle
		Vector3 dr = rf - ri;
		Vector3 pdr = center - ri;
		v = (normal*pdr)/(normal*dr); //intersection point
		if (v > 0.0 && v < 1.0 && !doub_equal(v*dr.magnitude(),0.0) ){return true;};
		return false;
		
  }

  //here r is a point on the surface
  __device__ __host__ Vector3 getNormal(Vector3 & r) const{
	return normal;
  }
		
  __device__ __host__ int getRegion(const Vector3 & r) const{
    return 0;
  }
		
  __device__ __host__  real getT2(const Vector3 & r) const{
	return -1.0;
  }
		
  __device__ __host__ real getD(const Vector3 & r) const{
	return -1.0;
  }
  
   __device__ __host__ real getPermeability() const{
	return 0.0;
  }
  
};
