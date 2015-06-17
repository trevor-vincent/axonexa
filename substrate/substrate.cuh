class Substrate {

private:

	real a, b, c;

public:

	Substrate(){}
	Substrate(real _a, real _b, real _c){
		a=_a; b=_b; c=_c;
	}

	__host__ __device__ real getA() const{return a;}
	__host__ __device__ real getB() const{return b;}
	__host__ __device__ real getC() const{return c;}

};
