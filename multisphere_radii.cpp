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
	
	//using namespace std;
	
	/* main routine */
	int main(int argc, char *argv[]){

		//cout.precision(30);

		/* 
		Variables From Command Line 
		Even if a certain geometry is selected, we will still declare the variables for the other geometries as well.
		So you will see variables for all geometries listed below, not just the one selected at the command line.
		*/

	FILELog::ReportingLevel() = logINFO;
		
	VARIABLES vars;
	vars = commandline_input(argc, argv);
	clock_t init, final;
		
	int number_of_particles = vars.num; 	
	int geometry = vars.geom;  // = 0 (Szafer Box Cells), = 1 (Square Packed Cylinders), = 2 (Hexagonally Packed Cylinders, = 3 (Randomly Packed Cylinders of random radii)
	double timestep = vars.dt; //Time interval (ms)
	double gspacing = vars.gs;
	bool use_T2 = vars.use_T2; // = false (no T2 decay), = true (T2 decay)


	int number_of_timesteps; number_of_timesteps = (int)ceil(gspacing/timestep); 
	int num_of_repeat = 4; //Number of times to repeat simulation. For every repeat all data is flushed and we start the simulation again.

	bool printparticleinfo = vars.printparticleinfo; //prints the trajectories of the particles if it is true
	bool printsignalinfo = vars.printsignalinfo; //prints the ADC calculations to the screen. This is useless, possibly now defunct variable.
	bool printinitialvars = vars.printinitialvars; //prints some interesting variables, but most defunct.


	//double t_Dz = (1-f)*D_extra + f*D_intra;
	//cout << "Theoretical Dz = " << t_Dz << endl;

	//SzaferBoxLattice lattice(D_extra, D_intra, T2_i, T2_e, cell_size_x, permeability, lattice_size);
		
	// Lattice<SquareCylinder_XY> lattice(D_extra, T2_e, permeability,1);
	// lattice.setLatticeVectors(lattice_size,lattice_size,lattice_size,Vector3(1.0,0.0,0.0),Vector3(0.0,1.0,0.0),Vector3(0.0,0.0,1.0));
	// lattice.addBasis(SquareCylinder_XY(lattice_size/2.0, lattice_size/2.0, cell_size_x, 1, T2_i, D_intra));

	//double width = cell_size_x;
	//Lattice<SquareCylinder_XY> lattice(D_extra, T2_e, permeability,1);
	//lattice.setLatticeVectors(lattice_size,lattice_size,lattice_size,Vector3(1.0,0.0,0.0),Vector3(0.0,1.0,0.0),Vector3(0.0,0.0,1.0));
	//lattice.addBasis(SquareCylinder_XY(lattice_size - width/2.0, lattice_size - width/2.0, width, 1, T2_i, D_intra));

	double f [] =  {0.6, 0.9,  0.7, 0.8};
	
	
	for (int kk = 0; kk < num_of_repeat; kk++) {
	
	double permeability = 0.0;
	double r_cytoplasm = .005;
	double r_nucleus = r_cytoplasm*f[kk];
	double D_extra = 1.82E-6;
	double D_nucleus = 1.31E-6 ;
	double D_cytoplasm = .48E-6;
	double T2 = 200;
	double a = .014996170545319518118759636223249;
	
	std::cout << "TRIAL " << kk << " NUCLEUS RADIUS = " << r_nucleus << std::endl;
		
	Lattice<MultiSphere> lattice(D_extra, T2, permeability);
	lattice.setLatticeVectors(a,a,a, Vector3(1.0,0.0,0.0), Vector3(0.0,1.0,0.0), Vector3(0.0,0.0,1.0));
	// lattice.addBasis(MultiSphere(0.5*a, 0.5*a, 0.5*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.0*a, 0.0*a, 0.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(1.0*a, 0.0*a, 0.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(1.0*a, 0.0*a, 1.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.0*a, 0.0*a, 1.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.0*a, 1.0*a, 0.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(1.0*a, 1.0*a, 0.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(1.0*a, 1.0*a, 1.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.0*a, 1.0*a, 1.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.5*a, 0.0*a, 0.5*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(1.0*a, 0.5*a, 0.5*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.5*a, 1.0*a, 0.5*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.0*a, 0.5*a, 0.5*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.5*a, 0.5*a, 1.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));
	lattice.addBasis(MultiSphere(0.5*a, 0.5*a, 0.0*a, r_nucleus, r_cytoplasm, T2,T2, D_nucleus, D_cytoplasm, 1 , 2));	
		
	int NOM = 10;
	int NOI = 100;
		
	Vector3 graddir1(1.0,0.0,0.0);
	//double graddir2 [] = {0.0, 1.0, 0.0};
	//double graddir3 [] = {0.0, 0.0, 1.0};
		
	OGSE_COS* measurements_x [NOM*NOI]; 

	double lnsignal [NOM];
	double b [NOM];
	

	
	
		for (int j = 0; j < NOI; j++){
		
			for (int i = 0; i < NOM; i++) {
					
				int N = 1+4*j;
				double G = i*0.0000025*N;
				double sigma = 40;	
				double gs = 2.0;
				double TE = 2.0*sigma + gs;
						
				measurements_x[i + j*NOM] = new OGSE_COS(N, sigma, gs, timestep, G, TE, number_of_particles, graddir1); 
					
			}
			
		}

		
		double ADCx, dt_eff, freq;

		
		
		Particles ensemble(number_of_particles,timestep, use_T2);
		lattice.initializeUniformly(ensemble.getGenerator(), ensemble.getEnsemble() );
		//lattice.initializeInBasis(ensemble.getGenerator(), ensemble.getEnsemble())
		//vector<int> regions = {2};
		//lattice.initializeInRegions(ensemble.getGenerator(), ensemble.getEnsemble(), regions)
		ensemble.printDistribution(0);
		//cout << " timestep = " << 0 << endl;
		
		for (int k = 0; k < NOM*NOI;k++){
		
			measurements_x[k]->updatePhase(ensemble.getEnsemble(), 0.0);
		
		}
		
		if (vars.printparticleposition == true) {ensemble.printparticlepositions(0.0, lattice);}
		
		for (int ii = 1; ii <= number_of_timesteps; ii++){
		
			if(ii % 10000 == 0) { 
			
				cout << "PERCENT DONE = "  << (ii*100)/number_of_timesteps << "%" << endl;	
				cout << "time = " << ii*timestep << ", Dx = " << ensemble.calcD(0, ii*timestep) << endl;
			
			}
				
			/* updates the position of all the particles */
			ensemble.updateposition(lattice);
			if (ii == number_of_timesteps) {ensemble.printDistribution(ii);}
			
			if (vars.printparticleposition == true) {ensemble.printparticlepositions(ii*timestep,lattice);}
			
			for (int k = 0; k < NOM*NOI; k++){
				
				measurements_x[k]->updatePhase(ensemble.getEnsemble(), ii*timestep);
				
			}  

		}
			
		
		for (int i = 0; i < NOI; i++){
		
			for (int j = 0; j < NOM; j++){
			freq = measurements_x[j + i*NOM]->get_freq();
			lnsignal[j] = log(measurements_x[j + i*NOM]->get_signal());
			b[j] = 	measurements_x[j + i*NOM]->get_b();
			//cout << lnsignal[j] << " " << b[j] << endl;
			}
			
			ADCx = linear_regression(lnsignal,b,NOM);
			ADCx = -1*ADCx; 
			cout << freq << " " << ADCx << endl;
		}
		
		for (int k = 0; k < NOM; k++){
		
			delete measurements_x[k];
			
		}  			

	}

		final= clock()-init; //final time - intial time
		cout << endl << "seconds for calculation: " << (double)final / ((double)CLOCKS_PER_SEC) << endl;
		cout << endl << " hours for calculation: " << (double)final/ (3600 * ((double)CLOCKS_PER_SEC)) << endl;

		if (vars.printtofile == true){ fclose(vars.fh);}
		if (vars.printparticleposition == true) {fclose(vars.fh_pp);}
		



		
		return 0;
	}




