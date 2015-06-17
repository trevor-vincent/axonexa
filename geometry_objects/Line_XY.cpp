#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include "log.h"
#include "compare_double.h"
#include "Objects.h"
#include "quadratic.h"	
#include <vector>
#include <algorithm>
#include <iostream>


Line_XY::Line_XY(double xi, double yi, double xf, double yf){

	_xi = xi;
	_yi = yi;

	_xf = xf;
	_yf = yf;

}


int Line_XY::getregion(Vector3 & r){
	return 0;
}

double Line_XY::getT2(Vector3 & r) {
	return -1.0;
}

double Line_XY::getD(Vector3 & r) {
	return -1.0;
}

Vector3 Line_XY::getNormal(double x, double y, double z){

	Vector3 zhat(0.0,0.0,1.0);
	Vector3 line( _xf - _xi, _yf - _yi, 0.0);
	
	Vector3 vprod = zhat % line; //vector product
	vprod.normalize();
	return vprod; 
	
}

bool Line_XY::intersect(Line_XY & ray2, double & v){
	
	double denom = ((_xf - _xi) * (ray2._yf - ray2._yi)) - ((_yf - _yi) * (ray2._xf - ray2._xi));
	double numer = ((_yi - ray2._yi) * (ray2._xf - ray2._xi)) - ((_xi - ray2._xi) * (ray2._yf - ray2._yi));
	double r = numer / denom;
	double numer2 = ((_yi - ray2._yi) * (_xf - _xi)) - ((_xi - ray2._xi) * (_yf - _yi));
	double s = numer2 / denom;

    if ((r < 0.0 || r > 1.0) || (s < 0.0 || s > 1.0)) {return false;}

	v = r;
	return true;

}

bool Line_XY::intersect(Line & ray2, double & v){
	
	double x1 = ray2.xi;
	double y1 = ray2.yi;
	double x2 = ray2.xf;
	double y2 = ray2.yf;
	
	double x3 = _xi;
	double y3 = _yi;
	double x4 = _xf;
	double y4 = _yf;
	
	double numerator = (x4-x3)*(y1-y3) - (y4 - y3)*(x1-x3);
	double denominator = (y4-y3)*(x2-x1) - (x4-x3)*(y2 - y1);
	
	if ( doub_equal(numerator, 0.0) || doub_equal(denominator,0.0) ){return false;}
	double r = numerator/denominator;
	double mag = ray2.magnitude();
    if (r > 0.0 && r < 1.0 && !doub_equal(r*mag, 0.0)) {v = r; return true;}
	return false;

}


