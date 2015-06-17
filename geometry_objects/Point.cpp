#include <ostream>
#include <limits>
#include <cmath>
#include "compare_double.h"
#include "Objects.h"
#include "quadratic.h"	
#include <vector>
#include <algorithm>
#include <iostream>



Point::Point(){

	x = 0.0;
	y = 0.0;
	z = 0.0;

}


Point::Point(double _x, double _y, double _z){

	x = _x;
	y = _y;
	z = _z;

}

Point::Point(double _x, double _y){

	x = _x;
	y = _y;
	z = 0.0;


}
Point::~Point(){
}
