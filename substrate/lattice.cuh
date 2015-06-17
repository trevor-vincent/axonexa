class Lattice {

private:

  real a, b, c;
  real T2_lattice;
  real T1_lattice;
  real D_lattice;
  int basisSize;

public:

  Lattice(){}
  Lattice(real _a, real _b, real _c, real _T2, real _T1, real _D, int _basisSize){
    T2_lattice = _T2;
	T1_lattice = _T1;
	D_lattice = _D;
	a=_a; b=_b; c=_c;
	basisSize = _basisSize;
  }

  __host__ __device__ real getA(){
	return a;
  }
  
  __host__ __device__ real getB(){
	return b;
  }
  
  __host__ __device__ real getC(){
	return c;
  }
  
  __device__ void initializeUniformly(Vector3 & r, curandState & localState) const
  {
    
	r.x = curand_uniform(&localState)*a;
	r.y = curand_uniform(&localState)*b;
	r.z = curand_uniform(&localState)*c;

  }
  
    __host__ void initializeUniformlyCPU(Vector3 & r) const
  {
    
    r.x = unifRandCPP()*a;
    r.y = unifRandCPP()*b;
    r.z = unifRandCPP()*c;

  }
  
  template <class Basis>
    __device__ void initializeInRegion( const Basis * basis, 
				     curandState & localState,
				     Vector3 & r,
					 int region
					 ) const

	  {
	  
		do {
		  r.x = curand_uniform(&localState)*a;
		  r.y = curand_uniform(&localState)*b;
		  r.z = curand_uniform(&localState)*c;
		} while ( inRegion(basis,r) != region);
	  
	  }
	  
  __device__ __host__ void correctBoundary(Vector3 & r) const
  {
    if (!inLatticeCell(r)){
      lmod(r);
    }
  }
  
  __device__ __host__ void lmod(Vector3 & r) const {

    if( r.x > a) { r.x = fmod(r.x,a); }
    if( r.x < 0.0) { r.x = fmod(r.x,a) + a; } 

    if( r.y > b) { r.y = fmod(r.y,b); }
    if( r.y < 0.0) { r.y = fmod(r.y,b) + b; } 

    if( r.z > c) { r.z = fmod(r.z,c); }
    if( r.z < 0.0) { r.z = fmod(r.z,c) + c; }       

  }
  
  template <class Basis>    
  __device__ __host__ double getT2(Basis * basis, Vector3 & r) const{

    double T2;
	Vector3 temp = r;
	lmod(temp);
    for (int i = 0; i < basisSize; i++){
      if( T2 = basis[i].getT2(temp), T2 > 0){
		return T2;
      }
    }
    return T2_lattice;

  }
  
   template <class Basis>
  __device__ __host__ double getT1(Basis * basis, Vector3 & r) const{

    double T1;
	Vector3 temp = r;
	lmod(temp);
    for (int i = 0; i < basisSize; i++){
      if( T1 = basis[i].getT1(temp), T1 > 0){
		return T1;
      }
    }
    return T1_lattice;

  }
  
 template <class Basis>
  __device__ __host__ double getD(const Basis * basis, const Vector3 & r) const{

    double D;
	Vector3 temp = r;
	lmod(temp);
    for (int i = 0; i < basisSize; i++){
      if( D = basis[i].getD(temp), D > 0){
	return D;
      }
    }

    return D_lattice;

  }
  
   template <class Basis>
  __device__ __host__ int inRegion(const Basis * basis, const Vector3 & r) const{

    int region;
	Vector3 temp = r;
	lmod(temp);
    for (int i = 0; i < basisSize; i++){
      if ( region = basis[i].getRegion(temp), basis[i].inside(temp)){
		return region;
      }
    }
    
    return 0;
  }

  
  __device__ __host__ bool inLatticeCell (const Vector3 & r) const{
    if ( r.x > a || r.x < 0.0 ||
			  r.y > b || r.y < 0.0 ||
					   r.z > c || r.z < 0.0 )
      {
	return false;
      }
    return true;
  }

   template <class Basis>
  __device__ __host__ bool intersectionCheckPermeable (const Basis *basis, 
					     const Vector3 & ri,
					     const Vector3 & rf,
					     real & v,
						 Vector3 & n,
						 real & permeability)const

  {

    //displacement vector

    Vector3 dr = rf - ri;
    
    //modded line 1

    Vector3 mod1_f = rf; lmod(mod1_f);
    Vector3 mod1_i = mod1_f - dr;
    
    //modded line 2
    Vector3 mod2_i = ri; lmod(mod2_i); 
    Vector3 mod2_f = mod2_i + dr;

    real vTemp = 0.0;
    real vBestEst = 10.0;

    for (int i = 0; i < basisSize; i++){

      if(basis[i].intersect(mod1_i, mod1_f, vTemp) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod1_i;
		n = basis[i].getNormal( intPoint );
		permeability = basis[i].getPermeability();
      }
      
      if(basis[i].intersect(mod2_i, mod2_f, vTemp) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod2_i;
		n = basis[i].getNormal( intPoint );
		permeability = basis[i].getPermeability();
      }

    } 

	v = vBestEst;
    return (vBestEst < 5.0);

  }
  
   template <class Basis>
  __device__ __host__ bool intersectionCheckImpermeable(const Basis *basis, 
					     const Vector3 & ri,
					     const Vector3 & rf,
					     real & v,
						 Vector3 & n)const

  {

    //displacement vector

    Vector3 dr = rf - ri;
    
    //modded line 1

    Vector3 mod1_f = rf; lmod(mod1_f);
    Vector3 mod1_i = mod1_f - dr;
    
    //modded line 2
    Vector3 mod2_i = ri; lmod(mod2_i); 
    Vector3 mod2_f = mod2_i + dr;

    real vTemp = 0.0;
    real vBestEst = 10.0;

    for (int i = 0; i < basisSize; i++){

      if(basis[i].intersect(mod1_i, mod1_f, vTemp) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod1_i;
		n = basis[i].getNormal( intPoint );
      }
      
      if(basis[i].intersect(mod2_i, mod2_f, vTemp) && vTemp < vBestEst){
		vBestEst = vTemp;
		Vector3 intPoint = dr*vBestEst + mod2_i;
		n = basis[i].getNormal( intPoint );
      }

    } 

	v = vBestEst;
    return (vBestEst < 5.0);

  }
  
  __device__ __host__ int getBasisSize(){
  
	return basisSize;
  
  }
  
  __device__ __host__ void setBasisSize(int _basisSize){
	basisSize = _basisSize;
  }

};
