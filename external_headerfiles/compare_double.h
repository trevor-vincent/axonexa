	/*  

	Name: compare_double.h
	Author: Trevor Vincent 
	Email: vincenttrevor@gmail.com
	Description:

	Floating point comparison (specificially double precision) comparison needs to be done by using subtraction and a tolerance
	These methods are used semi-frequently throughout the code and are thus made inline to speed up runtime.
	
	*/
	
	
	#ifndef _COMPARE_DOUBLE
	#define _COMPARE_DOUBLE

	//#define EPSILON std::numeric_limits<double>::epsilon()
	#define EPSILON 1e-14
	
	//is the double nan
	inline bool doub_isnan(double d)
	{
		/* standard is nan check */
		return d != d;
	}

	//smaller than or equal <=
	inline bool doub_stoe(double a, double b) {
		if (fabs(a - b) < EPSILON) {return true;}
		else if (a < b) {return true;}
		else {return false;}
	}

	inline bool doub_stoe(double a, double b, double tol) {
		if (fabs(a - b) < tol) {return true;}
		else if (a < b) {return true;}
		else {return false;}
	}

	//greater than or equal >=
	inline bool doub_gtoe(double a, double b) {
		if (fabs(a - b) < EPSILON) {return true;}
		else if (a > b) {return true;}
		else {return false;}
	}
	
	inline bool doub_gtoe(double a, double b, double tol) {
		if (fabs(a - b) < tol) {return true;}
		else if (a > b) {return true;}
		else {return false;}
	}

	// equals
	inline bool doub_equal(double a, double b) {
		if (fabs(a - b) < EPSILON) {return true;}
		else {return false;}
	}
	
	inline bool doub_equal(double a, double b, double tol) {
		if (fabs(a - b) < tol) {return true;}
		else {return false;}
	}


	#endif



