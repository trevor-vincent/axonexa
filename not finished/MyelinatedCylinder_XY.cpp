#include "Objects.h"
#include "quadratic.h"	
#include <vector>
#include <utility>
#include <iostream>

MyelinatedCylinder_XY::MyelinatedCylinder_XY(double center_x, double center_y, double radius_smallCylinder, double radius_bigCylinder, double T2_smallCylinder, double T2_bigCylinder, double D_smallCylinder, double D_bigCylinder, int region_smallCylinder, int region_bigCylinder){


	smallCylinder = Cylinder(center_x, center_y, radius_smallCylinder, T2_smallCylinder, D_smallCylinder, region_smallCylinder);
	bigCylinder = Cylinder(center_x, center_y, radius_bigCylinder,  T2_bigCylinder, D_bigCylinder, region_bigCylinder);

}


bool MyelinatedCylinder_XY::inside(double x []){

	return bigCylinder.inside(x);
	
}

int MyelinatedCylinder_XY::getregion(double x []){

	int region1 = smallCylinder.getregion(x);
	int region2 = bigCylinder.getregion(x);
	
	if ( region1 != 0 ){
		
		return region1;
		
	} 
	
	else if (region2 != 0){
	
		return region2;
	
	}

	return 0;

}

double MyelinatedCylinder_XY::getT2(double x []) {

	double T2_1 = smallCylinder.getT2(x);
	double T2_2 = bigCylinder.getT2(x);
	
	if ( T2_1 > 0.0){
		
		return T2_1;
		
	} 
	
	else if ( T2_2 > 0.0){
	
		return T2_2;
	
	}

	return -1.0;

}

double MyelinatedCylinder_XY::getD(double x []) {

	double D_1 = smallCylinder.getD(x);
	double D_2 = bigCylinder.getD(x);
	
	if ( D_1 > 0.0){
		
		return D_1;
		
	} 
	
	else if ( D_2 > 0.0){
	
		return D_2;
	
	}

	return -1.0;

}

bool MyelinatedCylinder_XY::intersect(Line & line, double & v){

	std::vector<double> pos_v;
	double v1 = 0.0, v2 = 0.0;
	
	if (smallCylinder.intersect(line, v1)){pos_v.push_back(v1);}
	if (bigCylinder.intersect(line, v2)){pos_v.push_back(v2);}
	if (pos_v.empty()){ return false;}
	
	std::sort(pos_v.begin(), pos_v.end());
	v = pos_v[0];
	return true;
}

Vector3 MyelinatedCylinder_XY::getNormal(double x, double y, double z){

	return smallCylinder.getNormal(x,y,z);

}

