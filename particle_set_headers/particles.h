	/*  

	Name: particles.h
	Author: Trevor Vincent 
	Email: vincenttrevor@gmail.com
	Description:

	The main algorithm used to update the particles position goes as:
	
	INITIALLY
	
	1) We uniformly distribute N particles over the lattice
	2) Determine the speeds for each of the particles
	3) Initialize all T2 degradation factors to 1.0
	
	DURING A TIMESTEP
	
	1) We generate a step vector with a random orientation and add this to the particles current position vector s
	2) We check if this step vector intersects any boundaries when added to the particles current position vector 
	3) Calculate T2 degradation factor 
	4) correct the position of the particle if it went outside the lattice
	5) update the speed of the particle if it went to a new region 
	
	REPEAT

	We use the Numerical Recipes Third Edition Random Generator Algorithm which has a period of 3.138x10^50 to make the random numbers

*/
	
	#ifndef PARTICLES_H_ //only allows one class declaration memory allotment to be made
	#define PARTICLES_H_

	class Particles  {

			
	
	private:
	
		double timestep; 
		bool use_T2;
		Ran randomnumber;
		vector<Particle> ensemble;


	public:

	//Constructors and destructor


	/***************************************************************************************
		What: Particles Constructor 
		Description:
		
		Initializes a Particles object, with the information sent from the command line and other areas.
		Also, a Ran (numerical recipes random number) struct is created and seeded with the current computer time
	 ****************************************************************************************/

	Particles(int num, double ts, bool useT2)  : randomnumber(time(0))  {

		use_T2 = useT2;
		timestep = ts;
		ensemble.resize(num);


	}
	
	Particles(int num, double ts, bool useT2, int seed) : randomnumber(seed)  {


		use_T2 = useT2;
		timestep = ts;
		ensemble.resize(num);


	}

	/*************************************************************************************** 
		What: Particles Destructor 
		Description:
		
		Kills a Particles object
	 ****************************************************************************************/

	~Particles(){}
	
	Ran& getGenerator(){
	
		return randomnumber;
	
	}
	
	
	vector<Particle>& getEnsemble(){
	
		return ensemble;
	
	}
	
	
	// template <class T>
	// void initialize_ensemble(Lattice<T> & lattice, int where){
	
		// //distribute particles uniformly everywhere in lattice
		// if(where == 0){lattice.initialize(randomnumber, ensemble);}
		
		// //distribute particles uniformly in Basis only
		// else if (where == 1){lattice.initializeInBasis(randomnumber, ensemble);}
		
		// //Soon there will be options for outside basis distribution and nonuniform distributions
		// else{std::cout << "OPTION NOT SUPPORTED YET, PARTICLES INITIALIZE ERROR" << std::endl; exit(0);}
		// updateSpeed(lattice);
	
	// }
	
	
	/***************************************************************************************
		What: updateposition method
		Description:
		
		updates the position of a point particle.
		
	 ****************************************************************************************/
	
	 
	 
	template <class T>
	void updateposition(Lattice<T> & lattice){
	
	updateSpeed(lattice);
	
	Vector3 last_position;
	Vector3 last_current_position;
	double theta,phi;
	double deltat;
	
	FILE_LOG(logDEBUG4) << endl << "****** UPDATE POSITION METHOD ********" << endl;

	for (int i=0; i < ensemble.size(); i++){

		deltat = 0.0;  
		
		last_position = ensemble[i].position;
		last_current_position = ensemble[i].position;
		
		FILE_LOG(logDEBUG4) << endl << " Before the Update " << endl;
		FILE_LOG(logDEBUG4) << " Particle " << i << endl << " " <<  ensemble[i] << endl; 
		FILE_LOG(logDEBUG4) << " Region = " << lattice.inregion(ensemble[i].position) << endl; 
		
		theta = acos(2*randomnumber.doub() - 1); //angle between 0 and 180 degrees
		phi = randomnumber.doub()*2*PI; //angle between 0 and 360 degrees
		
		ensemble[i].position += Vector3(ensemble[i].speed*timestep*sin(theta)*cos(phi),
										ensemble[i].speed*timestep*sin(theta)*sin(phi),
										ensemble[i].speed*timestep*cos(theta));

		FILE_LOG(logDEBUG4) << " After the Update " << endl;
		FILE_LOG(logDEBUG4) << " Particle " << i << endl << " " <<  ensemble[i] << endl; 
		FILE_LOG(logDEBUG4) << " Region = " << lattice.inregion(ensemble[i].position) << endl; 
		
		boundary(ensemble[i].position, last_current_position, last_position, deltat, lattice);
	
		FILE_LOG(logDEBUG4) << " After Boundary Call " << endl;
		FILE_LOG(logDEBUG4) << " Particle " << i << endl << " " <<  ensemble[i] << endl; 
		FILE_LOG(logDEBUG4) << " Region = " << lattice.inregion(ensemble[i].position) << endl; 
		
		ensemble[i].unbounded_position += (ensemble[i].position - last_position);
		
		FILE_LOG(logDEBUG4) << " CHANGED REGIONS? " << (lattice.inregion(ensemble[i].position)  != lattice.inregion(last_position)) << endl;
		// if (lattice.inregion(ensemble[i].position)  != lattice.inregion(last_position)){ 
			// std::cout << std::endl << "********PARTICLE CHANGED REGIONS*********" << std::endl;
			// std::cout << "Current Region = " << lattice.inregion(ensemble[i].position) << endl;
			// std::cout << "Last Region = " << lattice.inregion(last_position) << endl;
		// }
		
		if (use_T2 == true){

			double T2_i = lattice.getT2(last_position);
			double T2_f = lattice.getT2(ensemble[i].position);
			ensemble[i].T2factor *=exp(-deltat/T2_i)*exp((deltat - timestep)/T2_f);
			
		}

		
		lattice.correctBoundary(ensemble[i].position);
		
		FILE_LOG(logDEBUG4) << " AFTER CORRECT BOUNDARY " << endl;
		FILE_LOG(logDEBUG4) << " Particle " << i << endl << " " <<  ensemble[i] << endl; 
		FILE_LOG(logDEBUG4) << " Region = " << lattice.inregion(ensemble[i].position) << endl; 

		}
		
		/* After every step, the particle speeds are updated in case a particle has entered a new region */
		
	
		FILE_LOG(logDEBUG4) << endl << " ********************************************* " << endl;
	}


	/***************************************************************************************
		What: updateSpeed method
		Description:
		
		updates the speed of all the particles with a value that depends on their position.
		
		if, <r^2> = 6*D*t   (Eqn. 1)
		
		=> sqrt( < r^2 > ) = sqrt(6*D*t)   (Eqn. 2)
		
		Now divide both sides of (Eqn. 2) by t
		
		=> speed = sqrt(6*D/t)    

	 ****************************************************************************************/
	template <class T>
	void updateSpeed(Lattice<T> & lattice)
	{
	
		for (int i = 0; i < ensemble.size(); i++){
		
			ensemble[i].speed = sqrt(6*lattice.getD(ensemble[i].position)/timestep);
			
		}

	}

	/***************************************************************************************
		What: getdisplacement method
		Description:
		
		retrieves displacement along specified dimension dim
		dim = 0, x - axis
		dim = 1, y - axis
		dim = 2, z - axis
		
	 ****************************************************************************************/
 
	
	double calcD(int dim, double current_time){
	
		Vector3 difference;
		double meanxsqr = 0.0;


		
		for (int i = 0; i < ensemble.size(); i++){
			difference = ensemble[i].unbounded_position - ensemble[i].initial_position;
			if(dim ==0) {meanxsqr += difference.x*difference.x;}
			if(dim ==1) {meanxsqr += difference.y*difference.y;}
			if(dim ==2) {meanxsqr += difference.z*difference.z;}
		}

		meanxsqr /= ensemble.size();
		return (meanxsqr/(2.0*current_time));
	
	}
	
	template <class OBJECT>
	double calcDinRegion(int dim, double current_time, Lattice<OBJECT> lattice, int region){
	
		Vector3 difference;
		double meanxsqr = 0.0;
		
		for (int i = 0; i < ensemble.size(); i++){
		
			if (lattice.inregion(ensemble[i].position) == region){
				difference = ensemble[i].unbounded_position - ensemble[i].initial_position;
				if(dim == 0) {meanxsqr += difference.x*difference.x;}
				if(dim == 1) {meanxsqr += difference.y*difference.y;}
				if(dim == 2) {meanxsqr += difference.z*difference.z;}
			}
			
		}
		
		meanxsqr /= ensemble.size();
		return (meanxsqr/(2.0*current_time));
	
	
	}

	template <class T>
	void boundary(Vector3 & current_position, Vector3 & last_position, Vector3 & initial_position, double & deltat, Lattice<T> & lattice){//, bool reflection){

		FILE_LOG(logDEBUG4) << endl << "******* BOUNDARY METHOD *********" << endl;
		FILE_LOG(logDEBUG4) << "Current Position: " << current_position << endl; 
		FILE_LOG(logDEBUG4) << "Last Position: " << last_position << endl;
	
		bool transmitted = false;
		double current_speed = getSpeed(initial_position, lattice);
		double transmission_speed = getSpeed(current_position, lattice);
		double v = 0.0;
		Vector3 normal;
		
		while(lattice.intersection(current_position,last_position, normal, v) && transmitted == false){

			Vector3 dr = current_position - last_position;
			current_position = dr*v + last_position;
			deltat = deltat + (v*(dr.magnitude())/current_speed); 
			
			FILE_LOG(logDEBUG4) << endl << " AFTER INTERSECTION " << endl;
			FILE_LOG(logDEBUG4) << "Current Position: " << current_position << endl; 
			FILE_LOG(logDEBUG4) << "Last Position: " << last_position << endl;
			FILE_LOG(logDEBUG4) << "v: " << v << endl;
			FILE_LOG(logDEBUG4) << "deltat: " << deltat << endl;
		
			if (randomnumber.doub() > lattice.get_permeability()*4/current_speed) { 
				
				FILE_LOG(logDEBUG4) << endl << " REFLECTION" << endl;

				Vector3 incident = current_position - last_position;
				Vector3 reflected = incident - normal*2*(normal*incident);
				
				FILE_LOG(logDEBUG4) << "Normal: " << normal << endl;
				FILE_LOG(logDEBUG4) << "Incident: " << incident << endl; 
				FILE_LOG(logDEBUG4) << "Reflected: " << reflected << endl;
				
				reflected.normalize();
				last_position = current_position;
				current_position += reflected*current_speed*(timestep - (deltat));
				transmission_speed = getSpeed(current_position, lattice);
				
				FILE_LOG(logDEBUG4) << "Reflected Normalized: " << reflected << endl;
				FILE_LOG(logDEBUG4) << "Last Position " << last_position << endl;
				FILE_LOG(logDEBUG4) << "Current Position " << current_position << endl;
			}

			else { 
			
				Vector3 incident = current_position - last_position;
				incident.normalize();
				current_position += incident*transmission_speed*(timestep - (deltat));		
				transmitted = true;
			}

		} 
		
		FILE_LOG(logDEBUG4) << endl << " ********************************************* " << endl;
	}

	/***************************************************************************************
		What: getSpeed method #1
		Description:

		gets a particles speed based on its position
	 ****************************************************************************************/		
	template <class T>
	double getSpeed (Vector3 & position, Lattice <T> & lattice){

		return sqrt(6*lattice.getD(position)/timestep);
	}

	/***************************************************************************************
		What: printparticlepositions method
		Description:

		prints the particles trajectories to text files in the ./trajectories/ folder
		
		if unbounded = true, we print the unbounded position
		else, we print the bounded position
	****************************************************************************************/	
	template <class T>
	void printparticlepositions( double current_time ,  Lattice<T> & lattice) {


			for (int i = 0; i < ensemble.size(); i++ ) {
				
				stringstream num;
				string filename = "./trajectories/";
				FILE *fh;
				num << i; 
				filename = filename + num.str() + ".txt";
				
				if (doub_equal(current_time,0.0)) {

					fh = fopen((char*)filename.c_str(),"w");
					fprintf(fh, "%.14f %.14f %.14f %.14f %.14f %.14f %.5f %d %.14f \n", ensemble[i].position.x, ensemble[i].position.y,ensemble[i].position.z, ensemble[i].unbounded_position.x, ensemble[i].unbounded_position.y,ensemble[i].unbounded_position.z, current_time, lattice.inregion(ensemble[i].position), getSpeed(ensemble[i].position, lattice));		
					
				}
				
				else {
				
					fh = fopen((char*)filename.c_str(),"a");
					fprintf(fh, "%.14f %.14f %.14f %.14f %.14f %.14f %.5f %d %.14f \n", ensemble[i].position.x, ensemble[i].position.y,ensemble[i].position.z, ensemble[i].unbounded_position.x, ensemble[i].unbounded_position.y,ensemble[i].unbounded_position.z, current_time, lattice.inregion(ensemble[i].position), getSpeed(ensemble[i].position, lattice));		
				}
			
				fclose(fh);
				
			}
		
	
	}


	void printDistribution( int step ) {

		FILE *fh;
		stringstream num; num << step;
		string filename = "./distributions/" + num.str() + ".txt";


		for (int i = 0; i < ensemble.size(); i++ ) {
			
			
			if (i == 0){
				fh = fopen((char*)filename.c_str(),"w");
				fprintf(fh, "%.14f %.14f %.14f \n", ensemble[i].position.x, ensemble[i].position.y,ensemble[i].position.z);		
				fclose(fh);
			}
			
			else {
				fh = fopen((char*)filename.c_str(),"a");
				fprintf(fh, "%.14f %.14f %.14f \n", ensemble[i].position.x, ensemble[i].position.y,ensemble[i].position.z);		
				fclose(fh);			
			
			}
				
		}
		
	
	}	

	};
	#endif






