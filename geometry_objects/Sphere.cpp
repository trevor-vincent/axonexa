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


Sphere::Sphere(double _center_x, double _center_y, double _center_z, double _radius, double _T2, double _D, int _region){

	center_x = _center_x;
	center_y = _center_y;
	center_z = _center_z;
	radius = _radius;
	region = _region;
	T2 = _T2;
	D = _D;

}

Vector3 Sphere::getNormal(double x, double y, double z){

	double n_x = x - center_x; 
	double n_y = y - center_y;
	double n_z = z - center_z;
	double mag = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);
	
	
	// FILE_LOG(logDEBUG4) << "********SPHERE GETNORMAL METHOD********" << std::endl;
	// FILE_LOG(logDEBUG4) << "position = " << x << " " << y << " " << z << std::endl;
	// FILE_LOG(logDEBUG4) << "centers = " << center_x << " " << center_y << " " << center_z << std::endl;
	// FILE_LOG(logDEBUG4) << "n_x " << n_x << std::endl;
	// FILE_LOG(logDEBUG4) << "n_y " << n_y << std::endl;		
	// FILE_LOG(logDEBUG4) << "n_z " << n_z << std::endl;		
	
	
	return Vector3(n_x/mag, n_y/mag, n_z/mag);

}

bool Sphere::inside(Vector3 & r){

return (r.x - center_x)*(r.x - center_x) + (r.y - center_y)*(r.y - center_y) + (r.z - center_z)*(r.z - center_z) < radius*radius; 

}

int Sphere::getregion(Vector3 & r){

	if ( inside(r) ){
		
		return region;
		
	} 

	return 0;

}

double Sphere::getT2(Vector3 & r) {

	if ( inside(r) ){
		
		return T2;
		
	} 

	return -1.0;

}

double Sphere::getD(Vector3 & r) {

	if ( inside(r) ){
		
		return D;
		
	} 

	return -1.0;

}

bool Sphere::intersect(Line & line, double & v){

	std::vector<double> pos_v;
	double roots [2];
	double dx,dy,dz,l_x,l_y,l_z, a,b,c;
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
	l_z = center_z;
	a = dx*dx + dy*dy + dz*dz;
	b = 2.0*xi*dx - 2.0*dx*l_x + 2.0*yi*dy - 2.0*dy*l_y + 2.0*zi*dz - 2.0*dz*l_z;
	c = xi*xi + yi*yi +zi*zi -2*xi*l_x -2*yi*l_y -2*zi*l_z  + l_x*l_x + l_y*l_y + l_z*l_z - radius*radius;
	quadratic(a,b,c,roots);
	if(roots[0] > 0.0 && roots[0] < 1.0 && b*b>4*a*c && doub_equal(roots[0]*step_mag, 0.0) == false ){pos_v.push_back(roots[0]);}
	if(roots[1] > 0.0 && roots[1] < 1.0 && b*b>4*a*c && doub_equal(roots[1]*step_mag, 0.0) == false ){pos_v.push_back(roots[1]);}
	// FILE_LOG(logDEBUG4) << "********SPHERE INTERSECTION METHOD********" << std::endl;
	// FILE_LOG(logDEBUG4) << "center " << l_x << " " << l_y << " " << l_z << std::endl;
	// FILE_LOG(logDEBUG4) << "roots[0] " << roots[0] << std::endl;
	// FILE_LOG(logDEBUG4) << "roots[1] " << roots[1] << std::endl;		
	
	
	
	if (pos_v.empty()){ return false;}


	
	std::sort(pos_v.begin(), pos_v.end());
	v = pos_v[0];
	return true;
	
}



