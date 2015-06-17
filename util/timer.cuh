#include <time.h>

class GPUtimer {

	private:
	  bool mStarted, mStopped;
	  cudaEvent_t mStart, mStop;

	public:
	  GPUtimer() : mStarted(false), mStopped(false) {
	  	cudaEventCreate(&mStart);
		cudaEventCreate(&mStop);  
	  }
	  
	  ~GPUtimer() {
		cudaEventDestroy(mStart);
		cudaEventDestroy(mStop);
	  }
	  
	  void start(cudaStream_t s) { cudaEventRecord(mStart, s); 
									mStarted = true; mStopped = false; }
	  void stop(cudaStream_t s)  { assert(mStarted);
									   cudaEventRecord(mStop, s); 
									   mStarted = false; mStopped = true; }
	  float elapsed() {
		assert(mStopped);
		if (!mStopped) return 0; 
		cudaEventSynchronize(mStop);
		float elapsed = 0;
		cudaEventElapsedTime(&elapsed, mStart, mStop);
		return elapsed;
	  }
};

class CPUtimer {

	private:
	clock_t begin, end;
	double time_spent;
	
	public:
	CPUtimer(){};
	~CPUtimer(){};
	
	void start(){
		begin = clock();
	}

	void stop(){
		end = clock();
	}

	void display(){
		std::cout << std::setprecision(10);
		std::cout << " Elapsed CPU time = " << ((double) (end - begin)) / CLOCKS_PER_SEC << " s " << std::endl;
	}
	
	double getTime(){
		return ((double) (end - begin)) / CLOCKS_PER_SEC;
	}

};



