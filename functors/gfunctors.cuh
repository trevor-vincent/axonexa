class RectGrad{

	private:
	real gradient_strength;
	real gradient_duration;
	real gradient_spacing;
	real bvalue;
	Vector3 g_hat;
	
	public:
	__device__ __host__ RectGrad(){}
	__device__ __host__ RectGrad(real b, real gd, real gspac, Vector3 & ghat){
	
		bvalue = b;
		gradient_duration = gd;
		gradient_spacing = gspac;
		g_hat = ghat;
		gradient_strength = sqrt(b/(gradient_duration*gradient_duration*GAMMA*GAMMA*( (2.0/3.0)*gradient_duration + gradient_spacing)));
	
	}
	
	__device__ __host__ RectGrad(real gstr, real gd, real gspac){
	
		gradient_strength = gstr;
		gradient_duration = gd;
		gradient_spacing = gspac;
		g_hat = Vector3(1.0,0.0,0.0);
		bvalue = gradient_strength*gradient_strength*gradient_duration*gradient_duration*GAMMA*GAMMA*( (2.0/3.0)*gradient_duration + gradient_spacing);
	
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradient_duration) {
		
			return g_hat*gradient_strength;
		}
		
		else if ( time >= gradient_spacing + gradient_duration && time < 2.0*gradient_duration + gradient_spacing) {
		
			return g_hat*-gradient_strength;
		
		}
		
		else {
		
			return g_hat*0.0;
		}
	
	}
	
	__host__ real getB(){
	
		return bvalue;
	
	}
	
	__host__ real getDiffusionTime(){
	
		return gradient_spacing + gradient_duration;
	}
	
	__host__ real getG(){
	
		return gradient_strength;
	}
	
	__host__ real getTau(){
		
		return gradient_duration + gradient_spacing*.5;
	
	}

};


class CosGFunc{

	private:

  real G;
  real gradientDuration;
  real gradientSpacing;
  real period;
  real b;
	Vector3 gHat;
	
	public:
	__device__ __host__ CosGFunc(){}
	__device__ __host__ CosGFunc(real _G, real _gradientDuration, real _gradientSpacing, int _numOfPeriods, Vector3 & _gHat){
	
		
		gradientDuration = _gradientDuration;
		gradientSpacing = _gradientSpacing;
		gHat = _gHat;
		G = _G;
		period = gradientDuration/((real)_numOfPeriods);
		b = GAMMA*GAMMA*(1.0/4.0)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		else {
			return Vector3(0.0,0.0,0.0);
		}
	
	}
	
	//amplitude modifier
	__device__ __host__ real ampMod(real time) const{
	
		return cos(2.0*PI*time/period);
	}
	
	__device__ __host__ real getFreq(){return 1.0/period;}
	
	__device__ __host__ real getB(){
	
		return b;
	
	}
	
	__device__ __host__ real getG(){
	
		return G;
	
	}

};

class SinGFunc{

	private:

    real G;
    real gradientDuration;
    real gradientSpacing;
    real period;
    real b;
	Vector3 gHat;
	
	public:
	__device__ __host__ SinGFunc(){}
	__device__ __host__ SinGFunc(real _G, real _gradientDuration, real _gradientSpacing, int _numOfPeriods, Vector3 & _gHat){
	
		
		gradientDuration = _gradientDuration;
		gradientSpacing = _gradientSpacing;
		gHat = _gHat;
		G = _G;
		period = gradientDuration/((real)_numOfPeriods);
		b = GAMMA*GAMMA*(3.0/4.0)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		else {
			return Vector3(0.0,0.0,0.0);
		}
	
	}
	
	//amplitude modifier
	__device__ __host__ real ampMod(real time) const{
	
		return sin(2.0*PI*time/period);
	}
	
	__device__ __host__ real getFreq(){return 1.0/period;}
	
	__device__ __host__ real getB(){
		return b;
	}
	
	__device__ __host__ real getG(){
	
		return G;
	
	}

};

class ApodCosGFunc{

  private:

  real G;
  real gradientDuration;
  real gradientSpacing;
  real period;
  real b;
	Vector3 gHat;
	
	public:
	__device__ __host__ ApodCosGFunc(){}
	__device__ __host__ ApodCosGFunc(real _G, real _gradientDuration, real _gradientSpacing, int _numOfPeriods, Vector3 & _gHat){
		gradientDuration = _gradientDuration;
		gradientSpacing = _gradientSpacing;
		gHat = _gHat;
		G = _G;
		period = gradientDuration/((real)_numOfPeriods);
		b = GAMMA*GAMMA*(1.0/4.0)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods)*(1.0 - 1.0/(8.0*_numOfPeriods));
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		
		else {
			return Vector3(0.0,0.0,0.0);
		}
	
	}
	
	//amplitude modifier
	__device__ __host__ real ampMod(real time) const{
		if (time >= period/4.0 && time < gradientDuration - period/4.0 ){
			return cos(2.0*PI*time/period);
		} 

		else {
			return sin(4.0*PI*time/period);
		}	
	}
	
	__device__ __host__ real getFreq(){return 1.0/period;}
	__device__ __host__ real getB(){
		return b;
	}
	
	__device__ __host__ real getG(){
		return G;
	}

};


class SWOGSEFunc{

	private:

  real G;
  real gradientDuration;
  real gradientSpacing;
  real v; //frequency
  real b;
  Vector3 gHat;
	
	public:
	__device__ __host__ SWOGSEFunc(){}
	__device__ __host__ SWOGSEFunc(real _G, real _gradientDuration, real _gradientSpacing, int _numOfPeriods, Vector3 & _gHat){
	
		
		gradientDuration = _gradientDuration;
		gradientSpacing = _gradientSpacing;
		gHat = _gHat;
		G = _G;
		v = ((real)_numOfPeriods)/gradientDuration;
		double a1 = (1. - pow(-1., 2.*_numOfPeriods) + 4.*gradientDuration*v)/(4.*v);
		b = G*G*GAMMA*GAMMA*( gradientDuration/(6.*v*v) +  gradientSpacing*(gradientDuration - a1)*(gradientDuration - a1));
	}
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		else {
			return Vector3(0.0,0.0,0.0);
		}
	
	}
	
	// amplitude modifier
	__device__ __host__ real ampMod(real time) const{
		return pow(-1., 2.*time*v);
	}
	
	__device__ __host__ real getFreq(){return v;}
	
	__device__ __host__ real getB(){
	
		return b;
	
	}
	
	__device__ __host__ real getG(){
	
		return G;
	
	}

};

