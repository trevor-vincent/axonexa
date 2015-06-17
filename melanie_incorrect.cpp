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
	double D_extra = 2.5E-6;
	double D_intra = 1.0E-6;
	double T2_e = 200;
	double T2_i = 200;
		
	Vector3 xhat(1.0,0.0,0.0);
	Vector3 yhat(0.0,1.0,0.0);
	Vector3 zhat(0.0,0.0,1.0);		
		
	Lattice<SquareCylinder_XY> lattice(D_extra, T2_e, permeability);
	lattice.setLatticeVectors(lattice_size,lattice_size,lattice_size,xhat,yhat,zhat);
	lattice.addBasis(SquareCylinder_XY(lattice_size/2.0, lattice_size/2.0, width, 1, T2_i, D_intra));
	

	double ADCx,ADCy, ADCz, FA, mean_ADC, RA;
	
	double low_b = 108780;
	double big_b = 1000000;
	int num_of_dt = 8;
	double gspacings [] = { 8.0 , 10.5 , 14.0 , 18.5 , 24.5 , 32.5 , 42.5 , 56.5 };
	double G [2*num_of_dt];

	vector<double> low_b_sig_x(num_of_dt);
	vector<double> big_b_sig_x(num_of_dt);
	vector<double> low_b_sig_y(num_of_dt);
	vector<double> big_b_sig_y(num_of_dt);
	vector<double> low_b_sig_z(num_of_dt);
	vector<double> big_b_sig_z(num_of_dt);	
	
	vector<vector<double> > ADCx_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));
	vector<vector<double> > ADCy_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));
	vector<vector<double> > ADCz_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));
	vector<vector<double> > FA_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));
	vector<vector<double> > RA_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));
	vector<vector<double> > meanADC_rec(num_of_dt,vector<double>(num_of_repeat*num_of_dt));

	vector<vector<double> > ADCx_avg(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > ADCy_avg(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > ADCz_avg(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > FA_avg(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > RA_avg(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > meanADC_avg(num_of_dt,vector<double>(num_of_dt));
	
	vector<vector<double> > ADCx_err(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > ADCy_err(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > ADCz_err(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > FA_err(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > RA_err(num_of_dt,vector<double>(num_of_dt));
	vector<vector<double> > meanADC_err(num_of_dt,vector<double>(num_of_dt));

	for (int r = 0; r < num_of_repeat; r++) {

	vector<PGSE> measurements_x(2*num_of_dt);
	vector<PGSE> measurements_y(2*num_of_dt);
	vector<PGSE> measurements_z(2*num_of_dt);

	for (int i = 0; i < 8; i++){
		
		double echo_time = 2.0*grad_duration + gspacings[i];
		G[i] = sqrt(low_b/(GAMMA*GAMMA*grad_duration*grad_duration*(grad_duration + gspacings[i] - (grad_duration/3.0))));
		G[i + num_of_dt] = sqrt(big_b/(GAMMA*GAMMA*grad_duration*grad_duration*(grad_duration + gspacings[i] - (grad_duration/3.0))));
		
		measurements_x[i] = (PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, xhat));
		measurements_y[i] = (PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, yhat));
		measurements_z[i] = (PGSE(grad_duration,gspacings[i], timestep, G[0], echo_time, number_of_particles, zhat));
		
		measurements_x[i + num_of_dt] = (PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, xhat));
		measurements_y[i + num_of_dt] = (PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, yhat));
		measurements_z[i + num_of_dt] = (PGSE(grad_duration,gspacings[i], timestep, G[i+1], echo_time, number_of_particles, zhat));
	
	}
	
	vector<double> lnsignal(2);
	vector<double> b(2);
	
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
		
		
		for (int i = 0; i < num_of_dt; i++){
		
			fprintf(vars.fh_incorrectlinreg_ADCx, "\n low b time = %4.2f is constant \n", gspacings[i]);
			fprintf(vars.fh_incorrectlinreg_ADCy, "\n low b time = %4.2f is constant \n", gspacings[i]);
			fprintf(vars.fh_incorrectlinreg_ADCz, "\n low b time = %4.2f is constant \n", gspacings[i]);
			fprintf(vars.fh_incorrectlinreg_FA, "\n low b time = %4.2f is constant \n", gspacings[i]);
		
			for (int j = 0; j < num_of_dt; j++){
			double ADCx, ADCy, ADCz, RA, FA, mean_ADC, eff_diff_time;
			
			lnsignal[0] = log(measurements_x[i].get_signal());
			b[0] = 	measurements_x[i].get_b();
			lnsignal[1] = log(measurements_x[j+num_of_dt].get_signal());
			b[1] = 	measurements_x[j+num_of_dt].get_b();
			ADCx = -1.0*linear_regression(lnsignal,b);
			
			lnsignal[0] = log(measurements_y[i].get_signal());
			b[0] = 	measurements_y[i].get_b();
			lnsignal[1] = log(measurements_y[j+num_of_dt].get_signal());
			b[1] = 	measurements_y[j+num_of_dt].get_b();
			ADCy = -1.0*linear_regression(lnsignal,b);
			
			lnsignal[0] = log(measurements_z[i].get_signal());
			b[0] = 	measurements_z[i].get_b();
			lnsignal[1] = log(measurements_z[j+num_of_dt].get_signal());
			b[1] = 	measurements_z[j+num_of_dt].get_b();
			ADCz = -1.0*linear_regression(lnsignal,b);
			
			mean_ADC = (ADCx + ADCy + ADCz)/3;
			FA = sqrt(3.0/2.0)*sqrt((ADCx - mean_ADC)*(ADCx - mean_ADC) + (ADCy - mean_ADC)*(ADCy - mean_ADC) + (ADCz - mean_ADC)*(ADCz - mean_ADC))/sqrt(ADCx*ADCx + ADCy*ADCy + ADCz*ADCz);
			RA = sqrt(1.0/3.0)*sqrt((ADCx - mean_ADC)*(ADCx - mean_ADC) + (ADCy - mean_ADC)*(ADCy - mean_ADC) + (ADCz - mean_ADC)*(ADCz - mean_ADC))/mean_ADC;

			
			ADCx_rec[i][j + r*num_of_dt] = ADCx;
			ADCy_rec[i][j + r*num_of_dt] = ADCy;
			ADCz_rec[i][j + r*num_of_dt] = ADCz;
			FA_rec[i][j + r*num_of_dt] = FA;
			RA_rec[i][j + r*num_of_dt] = RA;
			meanADC_rec[i][j + r*num_of_dt] = mean_ADC;			
			
			
			eff_diff_time = measurements_x[i].get_DT();
			fprintf(vars.fh_incorrectlinreg_ADCx, "%4.2f %4.2f %.14f\n", gspacings[i], eff_diff_time, ADCx);
			fprintf(vars.fh_incorrectlinreg_ADCy, "%4.2f %4.2f %.14f\n", gspacings[i], eff_diff_time, ADCy);
			fprintf(vars.fh_incorrectlinreg_ADCz, "%4.2f %4.2f %.14f\n", gspacings[i], eff_diff_time, ADCz);
			fprintf(vars.fh_incorrectlinreg_FA, "%4.2f %4.2f %.14f %.14f %.14f\n", gspacings[i], eff_diff_time, FA, RA, mean_ADC);
			

			}
		}
			
	}
		
		fprintf(vars.fh_incorrectlinreg_ADCx, "\n MEAN & ERROR \n");
		fprintf(vars.fh_incorrectlinreg_ADCy, "\n MEAN & ERROR \n");
		fprintf(vars.fh_incorrectlinreg_ADCz, "\n MEAN & ERROR \n");
		fprintf(vars.fh_incorrectlinreg_FA, "\n MEAN & ERROR \n");
		for (int j = 0; j < num_of_dt; j++ ) {
		
			fprintf(vars.fh_incorrectlinreg_ADCx, "\n low b time = %4.2f is constant \n", gspacings[j]);
			fprintf(vars.fh_incorrectlinreg_ADCy, "\n low b time = %4.2f is constant \n", gspacings[j]);
			fprintf(vars.fh_incorrectlinreg_ADCz, "\n low b time = %4.2f is constant \n", gspacings[j]);
			fprintf(vars.fh_incorrectlinreg_FA, "\n low b time = %4.2f is constant \n", gspacings[j]);

			meanADC(ADCx_rec[j], num_of_dt, num_of_repeat, ADCx_avg[j], ADCx_err[j]);
			meanADC(ADCy_rec[j], num_of_dt, num_of_repeat, ADCy_avg[j], ADCy_err[j]);
			meanADC(ADCz_rec[j], num_of_dt, num_of_repeat, ADCz_avg[j], ADCz_err[j]);
			meanADC(FA_rec[j], num_of_dt, num_of_repeat, FA_avg[j], FA_err[j]);
			meanADC(RA_rec[j], num_of_dt, num_of_repeat, RA_avg[j], RA_err[j]);
			meanADC(meanADC_rec[j], num_of_dt, num_of_repeat, meanADC_avg[j], meanADC_err[j]);
			
			for (int i = 0; i < num_of_dt; i++ ) {

				fprintf(vars.fh_incorrectlinreg_ADCx, "%4.2f %.14f %.14f\n", gspacings[i], ADCx_avg[j][i], ADCx_err[j][i]);
				fprintf(vars.fh_incorrectlinreg_ADCy, "%4.2f %.14f %.14f\n", gspacings[i], ADCy_avg[j][i], ADCy_err[j][i]);
				fprintf(vars.fh_incorrectlinreg_ADCz, "%4.2f %.14f %.14f\n", gspacings[i], ADCz_avg[j][i], ADCz_err[j][i]);
				fprintf(vars.fh_incorrectlinreg_FA, "%4.2f %.14f %.14f %.14f %.14f %.14f %.14f\n", gspacings[i], FA_avg[j][i], FA_err[j][i], RA_avg[j][i], RA_err[j][i], meanADC_avg[j][i], meanADC_err[j][i]);
				
			}



		}
	

		final= clock()-init; //final time - intial time
		cout << endl << "seconds for calculation: " << (double)final / ((double)CLOCKS_PER_SEC) << endl;
		cout << endl << " hours for calculation: " << (double)final/ (3600 * ((double)CLOCKS_PER_SEC)) << endl;
		fclose(vars.fh_incorrectlinreg_ADCx);
		fclose(vars.fh_incorrectlinreg_ADCy);
		fclose(vars.fh_incorrectlinreg_ADCz);
		fclose(vars.fh_incorrectlinreg_FA);
		if (vars.printtofile == true){ fclose(vars.fh);}
		if (vars.printparticleposition == true) {fclose(vars.fh_pp);}	
		
		return 0;
	}




