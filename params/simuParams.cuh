class SimuParams {

public:

  int number_of_particles;
  int steps;
  int measurements;
  int seed;
  real timestep;
  real mx_initial;
  real my_initial;
  real mz_initial;

public:
  __host__ __device__ SimuParams(){}
  __host__ __device__ SimuParams(int _number_of_particles, 
				 int _steps, 
				 int _measurements, 
				 real _timestep, 
				 int _seed,
				 real _mx_initial,
				 real _my_initial,
				 real _mz_initial
				 )

  {
    number_of_particles = _number_of_particles;
    steps = _steps;
    measurements = _measurements;
    timestep = _timestep;
    seed = _seed;
    mx_initial = _mx_initial;
    my_initial = _my_initial;
    mz_initial = _mz_initial;


  }

};

class SimuParamsPhase {

public:

  int number_of_particles;
  int steps;
  int measurements;
  int seed;
  real timestep;
  real phase_initial;

public:
  __host__ __device__ SimuParamsPhase(){}
  __host__ __device__ SimuParamsPhase(
				int _number_of_particles, 
				 int _steps, 
				 int _measurements, 
				 real _timestep, 
				 int _seed,
				 real _phase_initial
				 )

  {
    number_of_particles = _number_of_particles;
    steps = _steps;
    measurements = _measurements;
    timestep = _timestep;
    seed = _seed;
    phase_initial = _phase_initial;
  }

};


class SimuParamsWC {

public:

  int number_of_particles;
  int steps;
  real timestep;
  real phase_initial;

public:
  __host__ __device__ SimuParamsWC(){}
  __host__ __device__ SimuParamsWC(
				 int _steps, 
				 real _timestep
				 )

  {
    steps = _steps;
    timestep = _timestep;
  }

};