	/***************************************************************************************
		What: Quadratic and Sgn Subroutines
		Author: Trevor Vincent (vincenttrevor@gmail.com)
		Description:
		
		The Quadratic Equation must be implemented in a certain way
		as to avoid great inaccuracy due to floating point arithmetic.
		
		For a more detailed explanation, read the corresponding section in Numerical Recipes or the 
		Wikipedia article on the quadratic equation.
		
		or the famous article:
		
		What Every Computer Scientist Should Know About Floating-Point Arithmetic, by David Goldberg,
		published in the March, 1991 issue of Computing Surveys. Copyright 1991,
		Association for Computing Machinery, Inc., reprinted by permission.
	 ****************************************************************************************/


	#ifndef _QUADRATIC
	#define _QUADRATIC

	/* returns the sign of a number (i.e is the number negative or positive) */
	inline double sgn(double x){
		if (x > 0) {return 1;}
		else if (x < 0) {return -1;}
		else {return 0;}
	}

	
	
	/* most accurate implementation of the quadratic equation */
	inline void quadratic(double a, double b, double c, double root[]){

		double q,x1,x2;

		q = -.5*(b + sgn(b)*sqrt(b*b - 4*a*c));
		x1 = q/a;
		x2 = c/q;

		root[0] = x1; //1st root
		root[1] = x2; //second root
		
	}

	#endif


