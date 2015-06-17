//Quadropolar Saddle Coil
class QSC_GRE{
	
public:
  real B0;
  real G;
  real gradientDuration;
  real gradientSpacing;
	
public:
  __device__ __host__ QSC_GRE(){
    B0 = 0.0;
    G = 0.0;
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
  }
	
  __device__ __host__ QSC_GRE(real _B0, real _G, real _gradientDuration, real _gradientSpacing){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
  }

  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
    if ( time < gradientDuration) {
      return Vector3( 0.0, G*r.z, B0 + G*r.y );
    }
		
    else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
      return Vector3( 0.0, -G*r.z, B0 - G*r.y );
    }
		
    else {
      return Vector3(0.0, 0.0, B0);
    }
	
  }
  
  __host__ real getB(){
 	return G*G*gradientDuration*gradientDuration*GAMMA*GAMMA*( (2.0/3.0)*gradientDuration + gradientSpacing);
  }
  
  __host__ real getDiffusionTime(){
	return gradientSpacing + gradientDuration;
  }
  
};


class QSC_PGSE{
	
public:
  real B0;
  real G;
  real gradientDuration;
  real gradientSpacing;
  real rfDuration;
  real dtmin;
  real dtmax;
  real TE;
  real B1;
	
public:
  __device__ __host__ QSC_PGSE(){
    B0 = 0.0;
    G = 0.0;
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
  }
	
  __device__ __host__ QSC_PGSE( real _B0, 
						   real _G, 
						   real _gradientDuration, 
						   real _gradientSpacing, 
						   real _rfDuration, 
						   real _dtmin, 
						   real _dtmax, 
						   real _TE
						  ){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	rfDuration = _rfDuration;
	dtmin = _dtmin;
	dtmax = _dtmax;
	TE = _TE;
	B1 = (PI)/fabs(GAMMA*rfDuration);
  }

  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
    if ( time < gradientDuration) {
      return Vector3(0.0, G*r.z, B0 + G*r.y);
    }

	else if ( time >= TE/2.0 && time < TE/2.0 + rfDuration ) {
	  real w = GAMMA*B0;
	  return Vector3( B1*cos(w*time), -B1*sin(w*time), B0);
	  // return Vector3( 0.0 , B1*cos(w*time), B0);
	}
	
    else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
      return Vector3(0.0, G*r.z, B0 + G*r.y);
    }
	
    else {
      return Vector3(0.0, 0.0, B0);
    }
	
  }
  
  __host__ real getB(){
 	return G*G*gradientDuration*gradientDuration*GAMMA*GAMMA*( (2.0/3.0)*gradientDuration + gradientSpacing);
  }
  
  __host__ real getDiffusionTime(){
	return gradientSpacing + gradientDuration;
  }
  
  __device__ real getStepTimeMin(real time){
     if (time >= TE/2.0 && time < TE/2.0 + rfDuration) {return dtmin;}
     else {return dtmax;}
  }
  
};


class QSC_COS_GRE{
	
public:
  real B0;
  real G;
  real gradientDuration;
  real gradientSpacing;
  real period;
  real b;  
  
public:
  __device__ __host__ QSC_COS_GRE(){
    B0 = 0.0;
    G = 0.0;
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
  }
	
  __device__ __host__ QSC_COS_GRE(real _B0, 
							real _G, 
							real _gradientDuration, real _gradientSpacing, real _numOfPeriods){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	b = GAMMA*GAMMA*(1.0/4.0)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);
	period = _gradientDuration/_numOfPeriods;
  }

  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
    if ( time < gradientDuration) {
      // return Vector3(0.0, 0.0, B0 + G*ampMod(time)*r.y);
	  return Vector3(0.0, G*ampMod(time)*r.z, B0 + G*ampMod(time)*r.y);

    }
		
    else if ( time >= gradientSpacing + gradientDuration && 
			  time < 2.0*gradientDuration + gradientSpacing) {
     
	return Vector3( 0.0, -G*ampMod(time - gradientSpacing - gradientDuration)*r.z, B0 - G*ampMod(time - gradientSpacing - gradientDuration)*r.y);
	 // return Vector3( 0.0, 0.0 , B0 - G*ampMod(time - gradientSpacing - gradientDuration)*r.y);
					
    }
		
    else {
      return Vector3(0.0, 0.0, B0);
    }
	
  }
  
  	
   __device__ __host__ real ampMod(real time) const{
	  return cos(2.0*PI*time/period);
   }

  __host__ real getB(){
 	return b;
  }
  
  __host__ real getFreq(){
	return 1.0/period;
  }
};



//with Trapezoidal pulses
class QSC_GRE_TRAP{
	
public:
  real B0;
  real G;
  real gradientDuration;
  real gradientSpacing;
  real gradientRampTime;
	
public:
  __device__ __host__ QSC_GRE_TRAP(){
    B0 = 0.0;
    G = 0.0;
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
  }
	
  __device__ __host__ QSC_GRE_TRAP(real _B0, real _G, real _gradientDuration, real _gradientSpacing, real _gradientRampTime){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	gradientRampTime = _gradientRampTime;
  }

  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
	real t1 = gradientRampTime;
	real t2 = gradientRampTime + gradientDuration;
	real t3 = 2.0*t1 + t2;
	real t4 = gradientSpacing;
	real t5 = t4 + t3;
	
	if (time < t1){
		real modif = ampMod1(time);
	   return Vector3(0.0, (G*r.z)*modif, B0 + (G*r.y)*modif);
	}
	
    else if ( time >= t1 && time < t2) {
      return Vector3(0.0, G*r.z, B0 + G*r.y);
    }
	
	else if ( time >= t2 && time < t3){
		real modif = ampMod2(time - t2);
		return Vector3(0.0, G*r.z*modif, B0 + G*r.y*modif);
	}

    else if ( time >= t5 && time < t5 + t1) {
		real modif = ampMod1(time - t5);
	   return Vector3(0.0, (G*r.z)*modif, B0 + (G*r.y)*modif);
    }
	
	else if ( time >= t5 + t1 && time < t5+t2 ){
      return Vector3(0.0, G*r.z, B0 + G*r.y);
	}
	
	else if ( time >= t5 + t2 && time < t5 + t3){
		real modif = ampMod2(time - t5 - t2);
		return Vector3(0.0, G*r.z*modif, B0 + G*r.y*modif);	
	}
	
    else {
      return Vector3(0.0, 0.0, B0);
    }
	
  }
  
  __device__ __host__ real ampMod1 (real time) const{
	return (G/gradientRampTime)*time;
  }
  
   __device__ __host__ real ampMod2 (real time) const{
	return G - (G/gradientRampTime)*time;
  }
  
  __host__ void printGz (real total_time, real timestep){
	int steps = (int)(timestep/total_time);
	for (int i = 0; i < steps; i++){
		std::cout << i*timestep << " " << operator()(Vector3(0.0,1.0,0.0), i*timestep) << std::endl;
	}
  }


};


class SIN_GRE{
	
public:
  real B0;
  Vector3 G;
  real gradientDuration;
  real gradientSpacing;
  real period;
	
public:
  __device__ __host__ SIN_GRE(){
    B0 = 0.0;
    G = Vector3(0.0,0.0,0.0);
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
	period = 0.0;
  }
	
	
  __device__ __host__ SIN_GRE(real _B0, Vector3 _G, real _gradientDuration, real _gradientSpacing, real _period){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	period = _period;
  }


  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
    if ( time < gradientDuration) {
      return Vector3(0.0, 0.0, B0 + G*r*sin(2.0*PI*time/period));
    }
		
    else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
      return Vector3(0.0, 0.0, B0 - G*r*sin(2.0*PI*(time - gradientSpacing - gradientDuration )/period) );
    }
		
    else {
      return Vector3(0.0, 0.0, B0);
    }
	
  }
	
};


class COS_GRE{
	
public:
  real B0;
  Vector3 G;
  real gradientDuration;
  real gradientSpacing;
  real period;
  real b;
	
public:
  __device__ __host__ COS_GRE(){
    B0 = 0.0;
    G = Vector3(0.0,0.0,0.0);
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
	period = 0.0;
  }
  __device__ __host__ COS_GRE(real _B0, Vector3 _G, real _gradientDuration, real _gradientSpacing, int _numOfPeriods){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	period = gradientDuration/((real)_numOfPeriods);
	b = GAMMA*GAMMA*(1.0/4.0)*G.magnitude()*G.magnitude()*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);

  }
    __device__ __host__ COS_GRE(real _B0, real _G, Vector3 G_hat, real _gradientDuration, real _gradientSpacing, int _numOfPeriods){
    B0 = _B0;
    G = G_hat*_G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
	period = gradientDuration/((real)_numOfPeriods);
	b = GAMMA*GAMMA*(1.0/4.0)*_G*_G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);

  }
  __device__ __host__ real getB(){
  
	return b;
	
  }
  __device__ __host__ real getFreq(){return 1.0/period;}
  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
    if ( time < gradientDuration) {
      return Vector3(0.0, 0.0, B0 + G*r*cos(2.0*PI*time/period));
    }
		
    else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
      return Vector3(0.0, 0.0, B0 - G*r*cos(2.0*PI*(time - gradientSpacing - gradientDuration )/period) );
    }	
    else {
      return Vector3(0.0, 0.0, B0);
    }
  }
	
};





class RECT_GRE{
public:
  real B0;
  Vector3 G;
  real gradientDuration;
  real gradientSpacing;
public:
  __device__ __host__ RECT_GRE(){
    B0 = 0.0;
    G = Vector3(0.0,0.0,0.0);
    gradientDuration = 0.0;
    gradientSpacing = 0.0;
  }
  __device__ __host__ RECT_GRE(real _B0, Vector3 _G, real _gradientDuration, real _gradientSpacing){
    B0 = _B0;
    G = _G;
    gradientDuration = _gradientDuration;
    gradientSpacing = _gradientSpacing;
  }
  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
    if ( time < gradientDuration) {
      return Vector3(0.0, 0.0, B0 + G*r);
    }
    else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
      return Vector3(0.0, 0.0, B0 - G*r);
    }	
    else {
      return Vector3(0.0, 0.0, B0);
    }
  }
};

class UCN_lowfield{
	
public:
  real B0;
  real G;
  real rfDuration;
  real dtmin;
  real dtmax;
  real w;
  real B1;
	
public:
  __device__ __host__ UCN_lowfield(){
    B0 = 0.0;
    G = 0.0;
  }
	
  __device__ __host__ UCN_lowfield( real _B0, 
						   real _G, 
						   real _rfDuration, 
						   real _dtmin, 
						   real _dtmax
						  ){
    B0 = _B0;
    G = _G;
	rfDuration = _rfDuration;
	dtmin = _dtmin;
	dtmax = _dtmax;
	B1 = (PI)/(2.0*fabs(GAMMA*rfDuration));
	w = GAMMA*B0;
  }

  __device__ __host__ Vector3 operator() (Vector3 r, real time) const {
	
	if ( time < rfDuration ) {
	  return Vector3( B1*cos(w*time) -.5*G*r.x, -B1*sin(w*time) -.5*G*r.y, G*r.z + B0);
	}
	
    else {
      return Vector3( -.5*G*r.x, -.5*G*r.y, G*r.z + B0);
    }
	
  }
  
  __device__ real getStepTimeMin(real time) const{
     if (time <= rfDuration) {return dtmin;}
     else {return dtmax;}
  }
  
};