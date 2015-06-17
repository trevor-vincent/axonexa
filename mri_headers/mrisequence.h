	#ifndef MRISEQUENCE_H_ 
	#define MRISEQUENCE_H_

	#define GAMMA 267522.2005
	#define PI 3.1415926535897932384626433832795

	class TannerStejkalSequence {
	
		protected:
		
		double b;
		double timestep;
		double grad_duration;
		double grad_strength;
		double grad_spacing;
		double echo_time;
		double signal;
		double diffusion_time;
	
		Vector3 graddir;
		std::vector<double> phase;
		
		public:
		
		TannerStejkalSequence(double _timestep, double _grad_duration, double _grad_strength, double _echo_time, double _grad_spacing, Vector3 _graddir, int number_of_particles){
		
			timestep = _timestep;
			grad_duration = _grad_duration;
			grad_strength = _grad_strength;
			grad_spacing = _grad_spacing;
			echo_time = _echo_time;
			graddir = _graddir;
			phase.resize(number_of_particles);
			std::fill( phase.begin(), phase.end(), 0.0 );
		
		}
		TannerStejkalSequence(){}
		
		double get_signal(){return signal;}
		double get_b(){return b;}
		double get_DT(){return diffusion_time;}
		
		void updatePhase(std::vector<Particle> & ensemble, double time){

			double gradmag = gradient(time);
		
			if (!doub_equal(gradmag,0.0) && !(time > echo_time) ){

				for (int i = 0; i < phase.size(); i++) {
				
					//phase[i] += GAMMA*gradmag*(graddir*ensemble[i].unbounded_position)*timestep;
					phase[i] += GAMMA*gradmag*(graddir*(ensemble[i].unbounded_position - ensemble[i].initial_position))*timestep;
				}
			
			}
			
			if ( doub_equal(time, echo_time, timestep/2.0) ) {
					
				calcSignal(ensemble);
					
			}

		}
		
		virtual double gradient(double time){
			if (doub_equal(time, 0.0, timestep/2.0)){
				return grad_strength;
			}
		
			else if (doub_equal(time, grad_spacing, timestep/2.0)){
				return -1.0*grad_strength;
			}
			
			else{
				return 0.0;
			}
		}
		
	
		void calcSignal(std::vector<Particle> & ensemble){
		
			double total_tmv_x = 0.0;
			double total_tmv_y = 0.0;
			
			//FILE_LOG(logDEBUG3) << std::endl << "***********UPDATE SIGNAL METHOD*********" << std::endl;	
			
			for (int i = 0; i < phase.size(); i++) {
			
				total_tmv_x += ensemble[i].T2factor*cos(phase[i]);
				total_tmv_y += ensemble[i].T2factor*sin(phase[i]);
					
				//FILE_LOG(logDEBUG3) << "ensemble[i].T2factor = " << ensemble[i].T2factor << std::endl;
				//FILE_LOG(logDEBUG3) << "total_tmv_x = " << total_tmv_x << std::endl;
				//FILE_LOG(logDEBUG3) << "total_tmv_y = " << total_tmv_y << std::endl;	

			}
			
			signal = (sqrt(total_tmv_x*total_tmv_x + total_tmv_y*total_tmv_y))/(double)phase.size();
		
		}
		
		//Signal in a boxed region with initial and final defining the diagonal of the box.
		void calcSignal(std::vector<Particle> & ensemble , Vector3 & initial, Vector3 & final){
		
			double total_tmv_x = 0.0;
			double total_tmv_y = 0.0;
			
			
			for (int i = 0; i < phase.size(); i++) {
			
				if( ensemble[i].unbounded_position > initial 
					&&
					ensemble[i].unbounded_position < final ){
					
					total_tmv_x += ensemble[i].T2factor*cos(phase[i]);
					total_tmv_y += ensemble[i].T2factor*sin(phase[i]);
					
				}
			}
			
			signal = (sqrt(total_tmv_x*total_tmv_x + total_tmv_y*total_tmv_y))/(double)phase.size();
		
		}
		
		template <class OBJECT>
		void calcSignal(std::vector<Particle> & ensemble, Lattice<OBJECT> & lattice, int region){
		
			double total_tmv_x = 0.0;
			double total_tmv_y = 0.0;
			
			//FILE_LOG(logDEBUG3) << std::endl << "***********UPDATE SIGNAL METHOD*********" << std::endl;	
			
			for (int i = 0; i < phase.size(); i++) {
			
				if (lattice.inregion(ensemble[i].position) == region){
					total_tmv_x += ensemble[i].T2factor*cos(phase[i]);
					total_tmv_y += ensemble[i].T2factor*sin(phase[i]);
				}
				
			}
			
			signal = (sqrt(total_tmv_x*total_tmv_x + total_tmv_y*total_tmv_y))/(double)phase.size();
		
		}
		
	
	};
	
	
	/**********************************************************************
	/*                             PGSE
	/*
	/*
	/*********************************************************************/	
	
	
	class PGSE : public TannerStejkalSequence {

		public:
		
		PGSE(double _grad_duration, double _grad_spacing, double _timestep, double _grad_strength, double _echo_time, int number_of_particles, Vector3 _graddir)
		:TannerStejkalSequence(_timestep,_grad_duration,_grad_strength,_echo_time,_grad_spacing, _graddir, number_of_particles)
		{
			diffusion_time = grad_duration + grad_spacing - (grad_duration/3.0);
			b = GAMMA*GAMMA*grad_strength*grad_strength*grad_duration*grad_duration*(grad_duration + grad_spacing - (grad_duration/3.0));
		}
		PGSE(){};
		~PGSE(){};

		double gradient(double t){
		
			double tol = timestep/2.0;

			
			if ( t < grad_duration ){

				return grad_strength;
			}
			
			else if ( doub_gtoe(t, grad_spacing + grad_duration, tol) && t < grad_spacing + 2.0*grad_duration) {
				return -grad_strength;
			}

			else{
			
				return 0.0;
			
			}
			
		}	


	
	};

	
	/**********************************************************************
	/*                             OGSE_SINE
	/*
	/*
	/*********************************************************************/
	class OGSE_SINE : public TannerStejkalSequence {

		private:


		int num_of_periods;
		double period;
		
		public:
		OGSE_SINE(int _num_of_periods, double _grad_duration, double _grad_spacing, double _timestep, double _grad_strength,  double _echo_time, int _number_of_particles, Vector3 _graddir)
		:TannerStejkalSequence(_timestep,_grad_duration,_grad_strength,_echo_time,_grad_spacing, _graddir, _number_of_particles){

			period = grad_duration/((double)_num_of_periods);
			b = GAMMA*GAMMA*(3.0/4.0)*grad_strength*grad_strength*grad_duration*grad_duration*grad_duration/(PI*PI*_num_of_periods*_num_of_periods);
			diffusion_time = (3.0/8.0)*period;
			
		}
		~OGSE_SINE(){};
		
		double gradient(double t){
			double tol = timestep/2.0;
			
			if ( t < grad_duration){
				return grad_strength*sin(2.0*PI*t/period);
			}
			
			else if ( doub_gtoe(t, grad_spacing + grad_duration, tol) && t < grad_spacing + 2.0*grad_duration) {
			  
				return -1.0*grad_strength*sin(2.0*PI*(t - grad_spacing - grad_duration )/period);
				
			}

			else{
				return 0.0;
			}
		}	
		
		double get_freq(){return 1.0/period;}
	
	};
	
	/**********************************************************************
	/*                             OGSE_COS
	/*
	/*
	/*********************************************************************/	
	class OGSE_COS  : public TannerStejkalSequence  {

		private:
		
		double period;

		public:
		OGSE_COS(int _num_of_periods, double _grad_duration, double _grad_spacing, double _timestep, double _grad_strength,  double _echo_time, int _number_of_particles, Vector3 _graddir)
		:TannerStejkalSequence(_timestep,_grad_duration,_grad_strength,_echo_time,_grad_spacing, _graddir, _number_of_particles){


			period = grad_duration/((double)_num_of_periods);
			b = GAMMA*GAMMA*(1.0/4.0)*grad_strength*grad_strength*grad_duration*grad_duration*grad_duration/(PI*PI*_num_of_periods*_num_of_periods);
			diffusion_time = (1.0/4.0)*period;
			
		}
		
		~OGSE_COS(){};

		double gradient(double t){
		
			double tol = timestep/2.0;
			
			if ( t < grad_duration){
				return grad_strength*cos(2.0*PI*t/period);
			}
			
			else if ( doub_gtoe(t, grad_spacing + grad_duration, tol) && t < grad_spacing + 2.0*grad_duration) {
			  
				return -1.0*grad_strength*cos(2.0*PI*(t - grad_spacing - grad_duration )/period);
				
			}

			else{
				return 0.0;
			}
		}	
		
		double get_freq(){return 1.0/period;}
		void get_details(std::ostream& out){
		
			out << std::endl << "OGSE COS SEQUENCE " << std::endl;
			out << " TE = " << echo_time << std::endl;
			out << " G = " << grad_strength << std::endl;
			out << " Sigma = " << grad_duration << std::endl;
			out << " Period = " << period << std::endl;
			out << " Frequency = " << (1.0/period) << std::endl;
			out << " Gradient Spacing = " << grad_spacing << std::endl;
			out << " Number of Periods = " << (grad_duration/period) << std::endl;
			out << " b = " << b << std::endl;
			out << " Graddir = " << graddir << std::endl;

		
		}
		
		friend std::ostream& operator<< (std::ostream &out, OGSE_COS & seq){
			seq.get_details(out);
			return out;
		}
	
	};	
	


	#endif
