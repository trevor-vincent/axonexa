/*

Evaluation of Bloch derivative without relaxation

*/

__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B){

  return (M%B)*GAMMA;

}

/*

Evaluation of Bloch derivative with relaxation

*/


__device__ __host__ Vector3 bloch(const Vector3 & M, const Vector3 & B, const real T1, const real T2){

  return (M%B)*GAMMA - Vector3(M.x/T2, M.y/T2, (M.z - 1.0f)/T1);

}

/*

Fourth Order Runge Kutta without Relaxation

*/


template <class FieldType>
__device__ __host__ void updateMagRK4(const Vector3 & r,
							   real* Mx, 
					           real* My, 
					           real* Mz, 
					           const FieldType & B, 
					           real i,
					           real h)
{

  real tn = i*h;
  Vector3 M(*Mx,*My,*Mz);
  Vector3 k1 = bloch(M, B(r,tn))*h;
  Vector3 k2 = bloch((M + k1*.5f), B(r,tn + .5f*h))*h;
  Vector3 k3 = bloch((M + k2*.5f), B(r,tn + .5f*h))*h;
  Vector3 k4 = bloch((M + k3),B(r,tn+h))*h;
  Vector3 finaldM = (k1 + k2*2.0f + k3*2.0f + k4)*(1.0f/6.0f);
  *Mx += finaldM.x;
  *My += finaldM.y;
  *Mz += finaldM.z;

}

/*

Fourth Order Runge Kutta withRelaxation

*/


template <class FieldType>
__device__ __host__ void updateMagRK4(real* Mx, 
					       real* My, 
					       real* Mz, 
					       const FieldType & B, 
						   real T1,
						   real T2,
					       real i,
					       real h)
{

  real tn = i*h;
  Vector3 M(*Mx,*My,*Mz);
  Vector3 k1 = bloch(M, B(tn), T1, T2)*h;
  Vector3 k2 = bloch((M + k1*.5f), B(tn + .5f*h), T1, T2)*h;
  Vector3 k3 = bloch((M + k2*.5f), B(tn + .5f*h), T1, T2)*h;
  Vector3 k4 = bloch((M + k3),B(tn+h), T1, T2)*h;
  Vector3 finaldM =(k1 + k2*2.0f + k3*2.0f + k4)*(1.0f/6.0f);
  *Mx += finaldM.x;
  *My += finaldM.y;
  *Mz += finaldM.z;

}
template <class FieldType>
__device__  void updateMagDORP(const Vector3 & M, const Vector3 & r, const FieldType & B, real tn, real & h, real hmax, real hmin, Vector3 & Mnew ){
    Vector3 k1 = bloch(M, B(r,tn))*h;
    Vector3 k2 = bloch((M + k1*.2f), B(r,tn + .2f*h))*h;
    Vector3 k3 = bloch((M + k1*(3.0f/40.0f) + k2*(9.0f/40.0f)), B(r,tn + .3*h))*h;
    Vector3 k4 = bloch((M + k1*(44.0f/45.0f) + k2*(-56.0f/15.0f) + k3*(32.0f/9.0f)),B(r,tn+.8*h))*h;
    Vector3 k5 = bloch((M + k1*(19372.0f/6561.0f) + k2*(-25360.0f/2187.0f) + k3*(64448.0f/6561.0f) + k4*(-212.0f/729)),B(r,tn+(8.0f/9.0f)*h))*h;
    Vector3 k6 = bloch((M + k1*(9017.0f/3168.0f) + k2*(-355.0f/33.0f) + k3*(-46732.0f/5247.0f) + k4*(49.0f/176.0f) + k5*(-5103.0f/18656.0f)),B(r,tn+h))*h;
    Vector3 k7 = bloch((M + k1*(35.0f/384.0f) + k3*(500.0f/1113.0f) + k4*(125.0f/192.0f) + k5*(-2187.0f/6784.0f) + k6*(11.0f/84.0f)),B(r,tn+h))*h;
	Mnew = M + k1*(35.0f/384.0f) + k3*(500.0f/1113.0f) + k4*(125.0f/192.0f) + k5*(-2187.0f/6784.0f) + k6*(11.0f/84.0f);
	Vector3 Mnew5 = M + k1*(5179.0f/57600.0f) + k3*(7571.0f/16695.0f) + k4*(393.0f/640.0f) + k5*(-92097.0f/339200.0f) + k6*(187.0f/2100.0f) + k7*(1.0f/40.0f);
	Vector3 diff = Mnew5-Mnew;
	Vector3 ratio; 
	ratio.x = diff.x/(.000001f + .000001f*max(Mnew5.x,Mnew.x));
    ratio.y = diff.y/(.000001f + .000001f*max(Mnew5.y,Mnew.y));
	ratio.z = diff.z/(.000001f + .000001f*max(Mnew5.z,Mnew.z));
	real err = sqrt(1.0f/3.0f*ratio.squareMagnitude());
    real hopt = h*(.98f)*pow( (1.0f/err), .2f);
	if ( hopt < hmax ){ h = hmin; }
}

template <class FieldType>
__device__ void updateMag_DORP(const Vector3 & r, real* Mx, real* My, real* Mz, const FieldType & B, int i, real h){

	real accumTime = 0.0f;
	real tk = i*h;
	real dt = h;
	real hmin = B.getStepTimeMin();

	Vector3 M(*Mx, *My, *Mz);
	Vector3 Mnew;
	
	do {

		updateMagDORP(M,r,B,tk,h,dt,hmin,Mnew);
		if ( h < (dt-accumTime) ){
		accumTime += h;
		updateMagDORP(M,r,B,tk,h,dt,hmin,Mnew);
		}
		
		else{
		accumTime += h;
		}
		
		tk = i*dt + accumTime;
		h = dt-accumTime;
		M = Mnew;
		
	} while ( accumTime < dt );

	*Mx = M.x;
	*My = M.y;
	*Mz = M.z;
	
}


__device__ __host__ Vector3 updateMag_noT2_rotate(Vector3 M, Vector3 B, real timestep){
 
  real Bmag = B.magnitude();
  Vector3 C = (B*(1.0f/Bmag))*(M*(B*(1.0f/Bmag)));
  Vector3 u = M-C;
  return (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0f/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
}

__device__ __host__ void updateMag_noT2_rotate(real* Mx, real* My, real* Mz, const Vector3 & B, real timestep){
  Vector3 M(*Mx, *My, *Mz);
  real Bmag = B.magnitude();
  Vector3 C = (B*(1.0f/Bmag))*(M*(B*(1.0f/Bmag)));
  Vector3 u = M-C;
  Vector3 finalM =  (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0f/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
  *Mx = finalM.x;
  *My = finalM.y;
  *Mz = finalM.z;
}
  

  
  
__host__ Vector3 updateMag_noT2_rotate_CPU_fast(Vector3 M, Vector3 B, real timestep){
  
  real Bmag = B.magnitude();
  Vector3 C = (B*(1.0f/Bmag))*(M*(B*(1.0f/Bmag)));
  Vector3 u = M-C;
  return (C + u*cos(-GAMMA*Bmag*timestep) + ((B*(1.0f/Bmag)) % u)*sin(-GAMMA*Bmag*timestep));
  
}


__device__ __host__ Vector3 updateMag_noT2_Euler(Vector3 M, Vector3 B, real timestep){
  return (M%B)*timestep*GAMMA;
}
