template <class T>
double bisec(T & functor, const double x1, const double x2, const double xacc){

	const int JMAX = 50;

	double dx, xmid, rtb;
	double f = functor(x1);
	double fmid = functor(x2);

	if (f*fmid >= 0.0) {

	std::cout << " BISECTION METHOD FAILED, ROOT ISN'T BRACKETED " << std::endl;
	exit(0);

	}

	rtb = f < 0.0 ? (dx = x2-x1, x1) : (dx=x1-x2,x2);
	for (int j =0; j < JMAX; j++){

		fmid = functor(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	
	std::cout << "ERROR: Too many bisections " << std::endl;
}