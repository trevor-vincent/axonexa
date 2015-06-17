//Line in the XY plane
class Line_XY {

public:

	double xi;
	double yi;

	double xf;
	double yf;
	
	Line_XY(double, double, double, double){};
	~Line_XY(){};
	bool intersect(Line_XY &, double &){};


}

Line_XY::Line_XY(double _xi, double _yi, double _xf, double _yf){

	xi = _xi;
	yi = _yi;

	xf = _xf;
	yf = _yf;

}

Line_XY::~Line_XY(){}


bool Line_XY::intersect(Line_XY & ray2, double & v){
	
	double denom = ((xf - xi) * (ray2.yf - ray2.yi)) - ((yf - yi) * (ray2.xf - ray2.xi));
	double numer = ((yi - ray2.yi) * (ray2.xf - ray2.xi)) - ((xi - ray2.xi) * (ray2.yf - ray2.yi));
	double r = numer / denom;
	double numer2 = ((yi - ray2.yi) * (xf - xi)) - ((xi - ray2.xi) * (yf - yi));
	double s = numer2 / denom;

    if ((r < 0.0 || r > 1.0) || (s < 0.0 || s > 1.0)) {return false;}

	v = r;
	return true;

}


