	/*  

	Name: test7.cpp
	Author: Trevor Vincent 
	Email: vincenttrevor@gmail.com
	Description:

	*/

	#include <iostream>
	#include <cmath>
	#include <sstream>
	#include <string>
	#include <limits>
	#include <vector>
	#include <algorithm>
	#include <utility>
	#include <string.h>
	#include <ctype.h>
	#include <ctime>
	#include <cstdlib>
	#include <limits>
	#include <cmath>
	
	#include "log.h"
	#include "compare_double.h"
	#include "getopt.h"	
	#include "Objects.h"
	#include "particle.h"
	#include "ran.h"
	#include "lattice.h"
	#include "particles.h"
	#include "mrisequence.h"
	#include "mean.h"
	#include "commandline_input.h"
	#include "bvalue.h"
	
	/* main routine */
	int main(int argc, char *argv[]){

	FILELog::ReportingLevel() = logINFO;
	VARIABLES vars;
	vars = commandline_input(argc, argv);
	clock_t init, final;
		
	int number_of_particles = vars.num; 
	double timestep = vars.dt;
	bool use_T2 = vars.use_T2; // = false (no T2 decay), = true (T2 decay)
	int num_of_repeat = vars.num_of_repeat; //Number of times to repeat simulation. For every repeat all data is flushed and we start the simulation again.
	int number_of_timesteps; number_of_timesteps = (int)ceil(vars.gs/timestep); 	
	
	double permeability = 0.0;
	double width = .00535;
	double lattice_size = 0.00601922026995422789215053328651;
	double grad_duration = 4.5;
	double D_extra = 2.0151414333359792147984537148051E-6;
	double D_intra = 1.0E-6;
	double T2_e = 200;
	double T2_i = 200;
		
	Vector3 xhat(1.0,0.0,0.0);
	Vector3 yhat(0.0,1.0,0.0);
	Vector3 zhat(0.0,0.0,1.0);		
		
	Lattice<SquareCylinder_XY> lattice(D_extra, T2_e, permeability);
	lattice.setLatticeVectors(lattice_size,lattice_size,lattice_size,xhat,yhat,zhat);
	//lattice.addBasis(SquareCylinder_XY(lattice_size/2.0, lattice_size/2.0, width, 1, T2_i, D_intra));
	


	double gspacings [] = { 8.0 , 10.5 , 14.0 , 18.5 , 24.5 , 32.5 , 42.5 , 56.5 };
	double bvals [] = {108780, 154720, 219040, 301730, 411980, 558990, 742420, 1000000 };
	double G [9];

	for (int kk = 0; kk < num_of_repeat; kk++) {

	vector<PGSE> measurements_x;
	vector<PGSE> measurements_y;
	vector<PGSE> measurements_z;

	for (int i = 0; i < 8; i++){
		double echo_time = 2.0*grad_duration + gspacings[i];
		G[0] = 0.0;
		G[i+1] = sqrt(bvals[i]/(GAMMA*GAMMA*grad_duration*grad_duration*(grad_duration + gspacings[i] - (grad_duration/3.0))));
		measurements_x.push_back(PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, xhat));
		measurements_y.push_back(PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, yhat));
		measurements_z.push_back(PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, zhat));
		measurements_x.push_back(PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, xhat));
		measurements_y.push_back(PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, yhat));
		measurements_z.push_back(PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, zhat));
	
	}
	
	vector<double> lnsignal(2);
	vector<double> b(2);
	

	
	cout << " trial = " << kk << endl;
		Particles ensemble(number_of_particles,timestep, use_T2);
		lattice.initializeUniformly(ensemble.getGenerator() , ensemble.getEnsemble() );
		
		for (int k = 0; k < measurements_x.size();k++){
			measurements_x[k].updatePhase(ensemble.getEnsemble(), 0.0);
			measurements_y[k].updatePhase(ensemble.getEnsemble(), 0.0);
			measurements_z[k].updatePhase(ensemble.getEnsemble(), 0.0);
		}
		
		for (int i = 1; i <= number_of_timesteps; i++){
			ensemble.updateposition(lattice);
			for (int k = 0; k < measurements_x.size();k++){
				measurements_x[k].updatePhase(ensemble.getEnsemble(), i*timestep);
				measurements_y[k].updatePhase(ensemble.getEnsemble(), i*timestep);
				measurements_z[k].updatePhase(ensemble.getEnsemble(), i*timestep);
			}
		}
		
		
		for (int i = 0; i < 8*2; i+=2){
			
			double ADCx, ADCy, ADCz, diff_time;
			
			lnsignal[0] = log(measurements_x[i].get_signal());
			b[0] = 	measurements_x[i].get_b();
			lnsignal[1] = log(measurements_x[i+1].get_signal());
			b[1] = 	measurements_x[i+1].get_b();
			ADCx = -1.0*linear_regression(lnsignal,b);
			
			std::cout << b[0] << " " << lnsignal[0] << " " << b[1] << " " << lnsignal[1] << " ";
			lnsignal[0] = log(measurements_y[i].get_signal());
			b[0] = 	measurements_y[i].get_b();
			lnsignal[1] = log(measurements_y[i+1].get_signal());
			b[1] = 	measurements_y[i+1].get_b();
			ADCy = -1.0*linear_regression(lnsignal,b);
			std::cout << b[0] << " " << lnsignal[0] << " " << b[1] << " " << lnsignal[1] << " ";
			lnsignal[0] = log(measurements_z[i].get_signal());
			b[0] = 	measurements_z[i].get_b();
			lnsignal[1] = log(measurements_z[i+1].get_signal());
			b[1] = 	measurements_z[i+1].get_b();
			ADCz = -1.0*linear_regression(lnsignal,b);
			std::cout << b[0] << " " << lnsignal[0] << " " << b[1] << " " << lnsignal[1] << std::endl;
			
			diff_time = measurements_x[i].get_DT();
			std::cout << i << " " << diff_time << " " << ADCx << " " << ADCy << " " << ADCz << " " << b[1] << " " << G[1] << std::endl; 
			
		}
			
	}
		

	

		final= clock()-init; //final time - intial time
		cout << endl << "seconds for calculation: " << (double)final / ((double)CLOCKS_PER_SEC) << endl;
		cout << endl << " hours for calculation: " << (double)final/ (3600 * ((double)CLOCKS_PER_SEC)) << endl;

		if (vars.printtofile == true){ fclose(vars.fh);}
		if (vars.printparticleposition == true) {fclose(vars.fh_pp);}	
		
		return 0;
	}




