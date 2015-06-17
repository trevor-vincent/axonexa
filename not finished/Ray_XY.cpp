class Ray_XY {

private:

	double xi;
	double yi;

	double xf;
	double yf;
	
	

public:

	Ray_XY(double, double, double, double){};
	~Ray_XY(){};
	bool rayXYrayXY_int(Ray_XY &, double &){};

}

Ray_XY::Ray_XY(double _xi, double _yi, double _xf, double _yf){

xi = _xi;
yi = _yi;

xf = _xf;
yf = _yf;

}

Ray_XY::~Ray_XY(){}


double Ray_XY::rayXYrayXY_int(Ray_XY & ray2, double & v){
	
	double denom = ((xf - xi) * (ray2.yf - ray2.yi)) - ((yf - yi) * (ray2.xf - ray2.xi));
	double numer = ((yi - ray2.yi) * (ray2.xf - ray2.xi)) - ((xi - ray2.xi) * (ray2.yf - ray2.yi));
	double r = numer / denom;
	double numer2 = ((yi - ray2.yi) * (xf - xi)) - ((xi - ray2.xi) * (yf - yi));
	double s = numer2 / denom;

    if ((r < 0.0 || r > 1.0) || (s < 0.0 || s > 1.0)) {return false;}

	v = r;
	return true;

}
