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
		gradient_strength = sqrt(b/(gradient_duration*gradient_duration*GAMMA*GAMMA*( (2.0f/3.0f)*gradient_duration + gradient_spacing)));
	
	}
	
	__device__ __host__ RectGrad(real gstr, real gd, real gspac){
	
		gradient_strength = gstr;
		gradient_duration = gd;
		gradient_spacing = gspac;
		g_hat = Vector3(1.0f,0.0f,0.0f);
		bvalue = gradient_strength*gradient_strength*gradient_duration*gradient_duration*GAMMA*GAMMA*( (2.0f/3.0f)*gradient_duration + gradient_spacing);
	
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradient_duration) {
		
			return g_hat*gradient_strength;
		}
		
		else if ( time >= gradient_spacing + gradient_duration && time < 2.0f*gradient_duration + gradient_spacing) {
		
			return g_hat*-gradient_strength;
		
		}
		
		else {
		
			return g_hat*0.0f;
		}
	
	}
	
	__host__ real getB(){
	
		return bvalue;
	
	}
	
	__host__ real getDiffusionTime(){
	
		return gradient_spacing + gradient_duration;
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
		b = GAMMA*GAMMA*(1.0f/4.0f)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0f*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		else {
			return Vector3(0.0f,0.0f,0.0f);
		}
	
	}
	
	//amplitude modifier
	__device__ __host__ real ampMod(real time) const{
	
		return cos(2.0f*PI*time/period);
	}
	
	__device__ __host__ real getFreq(){return 1.0f/period;}
	
	__device__ __host__ real getB(){
	
		return b;
	
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
		b = GAMMA*GAMMA*(1.0f/4.0f)*G*G*gradientDuration*gradientDuration*gradientDuration/(PI*PI*_numOfPeriods*_numOfPeriods);
	}
	
	__device__ __host__ Vector3 operator() (real time) const{
	
		if ( time < gradientDuration) {
			return (gHat*G)*ampMod(time);
		}
		
		else if ( time >= gradientSpacing + gradientDuration && time < 2.0f*gradientDuration + gradientSpacing) {
			return (gHat*-G)*ampMod(time - gradientSpacing - gradientDuration);
		}
		else {
			return Vector3(0.0f,0.0f,0.0f);
		}
	
	}
	
	//amplitude modifier
	__device__ __host__ real ampMod(real time) const{
	
		return sin(2.0f*PI*time/period);
	}
	
	__device__ __host__ real getFreq(){return 1.0f/period;}
	
	__device__ __host__ real getB(){
	
		return b;
	
	}

};