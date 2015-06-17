#include <ostream>
#include <limits>
#include <cmath>
#include "compare_double.h"
#include "Objects.h"
#include "quadratic.h"	
#include <vector>
#include <algorithm>
#include <iostream>


//3d line
Line::Line(double _xi, double _yi, double _zi, double _xf, double _yf, double _zf){

	xi = _xi;
	yi = _yi;
	zi = _zi;

	xf = _xf;
	yf = _yf;
	zf = _zf;


}

//2d line on the XY plane
Line::Line(double _xi, double _yi, double _xf, double _yf){

	xi = _xi;
	yi = _yi;

	xf = _xf;
	yf = _yf;
	
	zi = 0.0;
	zf = 0.0;

}



Line::~Line(){};

void Line::operator *= (const double v){
		
		xf = (xf - xi)*v + xi;
		yf = (yf - yi)*v + yi;
		zf = (zf - zi)*v + zi;
			
}

double Line::magnitude(){

	return sqrt( (xf-xi)*(xf-xi) + (yf-yi)*(yf-yi) + (zf-zi)*(zf-zi) );

}

Line Line::operator*(const double v) const{
			
	return Line(xi,yi,zi,(xf - xi)*v + xi,(yf - yi)*v + yi,(zf - zi)*v + zi);
}
