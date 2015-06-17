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



Cylinder_XY::Cylinder_XY(double _center_x, double _center_y, double _radius, double _T2, double _D, int _region){

	center_x = _center_x;
	center_y = _center_y;
	radius = _radius;
	region = _region;
	T2 = _T2;
	D = _D;
	
}


bool Cylinder_XY::inside(Vector3 & r){

	return ((r.x - center_x)*(r.x - center_x) + (r.y - center_y)*(r.y - center_y) < radius*radius);

}

int Cylinder_XY::getregion(Vector3 & r){

	if ( inside(r) ){
		
		return region;
		
	} 

	return 0;

}

double Cylinder_XY::getT2(Vector3 & r) {

	if ( inside(r) ){
		
		return T2;
		
	} 

	return -1.0;

}

double Cylinder_XY::getD(Vector3 & r) {

	if ( inside(r) ){
		
		return D;
		
	} 

	return -1.0;

}

Vector3 Cylinder_XY::getNormal(double x, double y, double z){

	double n_x = x - center_x; 
	double n_y = y - center_y;
	double mag = sqrt(n_x*n_x + n_y*n_y);
	return Vector3(n_x/mag, n_y/mag, 0.0);

}


bool Cylinder_XY::intersect(Line & line, double & v){

	std::vector<double> pos_v;
	double roots [2];
	double dx,dy,dz,l_x,l_y, a,b,c;
	double xi = line.xi;
	double xf = line.xf;
	double yi = line.yi;
	double yf = line.yf;
	double zi = line.zi;
	double zf = line.zf;
	
	dx = xf - xi;
	dy = yf - yi;
	dz = zf - zi;

	double step_mag = sqrt(dx*dx + dy*dy + dz*dz);	

	l_x = center_x;
	l_y = center_y;
	a = dx*dx + dy*dy;
	b = 2*xi*dx - 2*dx*l_x + 2*yi*dy - 2*dy*l_y;
	c = xi*xi + yi*yi -2*xi*l_x -2*yi*l_y + l_x*l_x + l_y*l_y - radius*radius;
	quadratic(a,b,c,roots);
	
	if(roots[0] > 0.0 && roots[0] < 1.0 && b*b>4*a*c && doub_equal(roots[0]*step_mag, 0.0) == false ){pos_v.push_back(roots[0]);}
	if(roots[1] > 0.0 && roots[1] < 1.0 && b*b>4*a*c && doub_equal(roots[1]*step_mag, 0.0) == false ){pos_v.push_back(roots[1]);}

	
	if (pos_v.empty()){ return false;}

	
	std::sort(pos_v.begin(), pos_v.end());
	v = pos_v[0];
	
	
	return true;
	
}



