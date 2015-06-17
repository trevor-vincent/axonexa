class iGaussFunc {

	private:
	double _C;
	
	public:
	diffFunc(double C){	_C = C;}
	~diffFunc(){};

	double operator() (double R){
	
		double pie = 3.1415926535897932384626433832795;
		return ( erf(R) - R*exp(-1.0*R*R)*2.0/(sqrt(pie)) - _C);
	
	
	}
};

class stepGenerator{
m

	private:
	
		static int _option;
		static int _seed;
		static Ran _randomnumber;
		
		static vector<Vector3> _dir;
		static vector<double> _gauss;

	
	public:
	stepGenerator(){}
	~stepGenerator(){}	
	
	static void setMembers(int seed, int option){
		_option = option;
		_seed = seed;
		_randomnumber = Ran(seed);
	}
	
	static void createGaussTable(int length){
		
		double inc = 1.0/length;
		_gauss.resize(length);
		_gauss[0] = 0.0;
		
		for (int i = 1; i < length; i++){
			
			double c = i*inc;
			double x1 = _gauss[i-1];
			double x2 = (x1+1.0)*2.0;
			double xacc = (1E-15)*(fabs(x1) + fabs(x2));
			iGaussFunc functor(c);
			_gauss[i] = bisec(functor, x1, x2, xacc);
			
		}
		
	}
	
	//generates vectors in only an octant of the sphere using Marsaglia's method
	static void createVectorTable(int length){
	
		double x1,x2,x,y,z;
		_dir.resize(length);
		for (int i = 0; i < length; i++){
		
			x1 = 2.0*(_randomnumber.doub()) - 1.0; 
			x2 = 2.0*(_randomnumber.doub()) - 1.0;
			
			if ( x1*x1 + x2*x2 >= 1.0) {
				
				x = 2*x1*sqrt(1-x1*x1-x2*x2);
				y = 2*x2*sqrt(1-x1*x2-x2*x2);
				z = 1-2(x1*x1 + x2*x2);
				
				if ( x > 0 && y > 0 && z > 0){_dir[i] = Vector3(x,y,z);}
				else {i--;}
			
			}
			
			else {i--;}
		
		}
	
	}
		
	static double generateSpeed(double D, double dt){
		
		
		//Gaussian table
		if (option == 0){

			double length = _gauss.size();
			double inc = 1.0/length;
			double X = _randomnumber.doub();
			int i = floor( X / inc );
			
			return R*sqrt(6*D/dt);

		}
		
		//Fixed velocity ala Szafer & Camino
		else {
		
			return sqrt(6*D/dt);
		
		}
		
		
	}
	
	//D = diffusion coeff, v = speed, dt = timestep
	static Vector3 generateStep(double D, double v, double dt){
	
		//Gaussian table
		if (option == 0){

			int dev = (_randomnumber.int64() % (_dir.size() - 1) );
			if(_randomnumber.doub() > .5){_dir[dev].x*(-1.0);}
			if(_randomnumber.doub() > .5){_dir[dev].y*(-1.0);}
			if(_randomnumber.doub() > .5){_dir[dev].z*(-1.0);}
			return (_dir[dev]*v); 
			
		}
		
		//Fixed step ala Szafer & Camino
		else {
		
			return Vector3(v*timestep*sin(theta)*cos(phi),
						   v*timestep*sin(theta)*sin(phi),
						   v*timestep*cos(theta));		
			
		}		
	
	}

		

};

