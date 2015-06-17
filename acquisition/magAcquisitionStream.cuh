//add devstates to kernel stream
template <class FieldType>
class magAcquisitionStream {

private:
	
  int steps_max; 
  int meas_max;

  //put this into simuparams and remove it from here
  int number_of_particles;
		
public:
	
  std::vector<magAcquisition<FieldType> > acqSet;

  magAcquisitionStream(int num){number_of_particles = num; steps_max = -1; meas_max = -1;}
  ~magAcquisitionStream(){}
		
  void addAcquisition(magAcquisition<FieldType> acq){
    acqSet.push_back(acq);
    if (acq.getNumOfSteps() > steps_max){ steps_max = acq.getNumOfSteps(); }
    if (acq.getNumOfMeas() > meas_max){ meas_max = acq.getNumOfMeas(); }
  }
  
  magAcquisition<FieldType> & getAcquisition(int num){
  
    return acqSet[num];
  
  }
		
  //put in Acquisition class
  template <class Basis>
  void runAcquisition(int num, Basis & basis, Vector3 & initialM, real timestep, int blocks, int threads, int devNum, int numOfSM){
	
    int steps = acqSet[num].getNumOfSteps();
    int measurements = acqSet[num].getNumOfMeas();

	
    cudaVector<FieldType> dev_fields(acqSet[num].getFieldFunctors().size());
    cudaVector<real> dev_mx(measurements*number_of_particles);
    cudaVector<real> dev_my(measurements*number_of_particles);
    cudaVector<real> dev_mz(measurements*number_of_particles);
    
	cudaVector<curandState> dev_states(number_of_particles);
	
	cudaVector<real> dev_total_mx(measurements);
	cudaVector<real> dev_total_my(measurements);
	cudaVector<real> dev_total_mz(measurements);
	
    cudaScalar<Basis> dev_basis( basis );
    cudaScalar<SimuParams> dev_par( SimuParams(number_of_particles,
					       steps,
					       measurements, 
					       timestep, 
					       acqSet[num].getSeed(),
					       initialM.x, 
					       initialM.y, 
					       initialM.z ) ); 

    dev_fields = acqSet[num].getFieldFunctors();
    dev_fields.copyToDevice();
    dev_basis.copyToDevice();
    dev_par.copyToDevice();

	setup_kernel <<< blocks, threads >>> (dev_states.getPointer(), dev_par.getPointer());
	
	 updateWalkersMag<FieldType, Basis> <<< blocks, threads >>> (

								    dev_par.getPointer(), 
								    dev_basis.getPointer(),  
								    dev_fields.getPointer(), 
								    dev_states.getPointer(),
								    dev_mx.getPointer(), 
								    dev_my.getPointer(), 
								    dev_mz.getPointer()
								      
								    );						  

	// cudaLastErrorCheck();
	// std::cout << "ERROR HERE?" << std::endl;
									
	dev_mx.sum(dev_total_mx, 1024, numOfSM, measurements, number_of_particles, 0);
	dev_my.sum(dev_total_my, 1024, numOfSM, measurements, number_of_particles, 0);
	dev_mz.sum(dev_total_mz, 1024, numOfSM, measurements, number_of_particles, 0);

	dev_total_mx.copyFromDevice();
    dev_total_my.copyFromDevice();
    dev_total_mz.copyFromDevice();
	
    dev_total_mx.copyTo(acqSet[num].getMx());
    dev_total_my.copyTo(acqSet[num].getMy());
    dev_total_mz.copyTo(acqSet[num].getMz());

  }

    template <class Basis>
    void runCPUAcquisition(int num, Basis & basis, Vector3 & initialM, real timestep){
  
	
    int steps = acqSet[num].getNumOfSteps();
    int measurements = acqSet[num].getNumOfMeas();
 
	std::vector<real> mx(number_of_particles*measurements);
	std::vector<real> my(number_of_particles*measurements);
	std::vector<real> mz(number_of_particles*measurements);
	
	
	
	std::vector<real> total_mx(measurements);
	std::vector<real> total_my(measurements);
	std::vector<real> total_mz(measurements);
 
	SimuParams pars(number_of_particles,
				steps,
				measurements, 
				timestep, 
				acqSet[num].getSeed(),
				initialM.x, 
				initialM.y, 
				initialM.z ); 	

	srand(acqSet[num].getSeed());
	
	updateWalkersMagCPU<FieldType, Basis>(&pars, &basis, &acqSet[num].getFieldFunctors()[0], &mx[0], &my[0], &mz[0]);
	
	for (int i = 0; i < measurements; i++){
		
		real total_mx_temp = 0.0;
		real total_my_temp = 0.0;
		real total_mz_temp = 0.0;
		
		for (int j = 0; j < number_of_particles; j++){
		
			total_mx_temp += mx[j + i*number_of_particles];
			total_my_temp += my[j + i*number_of_particles];
			total_mz_temp += mz[j + i*number_of_particles];
			
		}
	
		total_mx[i] = total_mx_temp;
		total_my[i] = total_my_temp;
		total_mz[i] = total_mz_temp;
	
	}
	
	acqSet[num].getMx() = total_mx;
	acqSet[num].getMy() = total_my;
	acqSet[num].getMz() = total_mz;
	
  } 
  
template <class Basis>		
  void runAcquisitionStream(Basis & basis, 
							Vector3 & initialM, 
							real timestep, 
							int blocks, 
							int threads , 
							int numOfDevices, 
							std::vector<int> & plan,
							std::vector<int> & numOfSMPerDevice){
	
    vector< pinnedVector<FieldType> > host_fields(acqSet.size()); 
	vector< pinnedScalar<SimuParams > > host_pars(acqSet.size());
    vector< pinnedVector<real> > host_mx(acqSet.size()); 
	vector< pinnedVector<real> > host_my(acqSet.size()); 
	vector< pinnedVector<real> > host_mz(acqSet.size()); 

    for (int i = 0; i < acqSet.size(); i++){
      
		host_fields[i].alloc(meas_max);
		host_mx[i].alloc(meas_max);
		host_my[i].alloc(meas_max);
		host_mz[i].alloc(meas_max);

    } 
	
	
    std::vector< std::vector< cudaVector< FieldType > > > dev_fields(numOfDevices,std::vector< cudaVector< FieldType > >(2) );
	std::vector< std::vector< cudaScalar<SimuParams> > > dev_pars(numOfDevices, std::vector< cudaScalar< SimuParams > >(2) );
	std::vector< std::vector< cudaVector<real> > >dev_mx(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<real> > >dev_my(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<real> > >dev_mz(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<real> > >dev_total_mx(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<real> > >dev_total_my(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<real> > >dev_total_mz(numOfDevices,std::vector< cudaVector<real> >(2));
	std::vector< std::vector< cudaVector<curandState> > > devStates(numOfDevices,std::vector< cudaVector<curandState> >(2));
	std::vector< std::vector< cudaScalar<Basis> > > dev_basis(numOfDevices, std::vector< cudaScalar< Basis > >(2) );
	std::vector < std::vector<cudaStream_t> >streams(numOfDevices, std::vector<cudaStream_t>(2));

	// CPUtimer ctimer; ctimer.start();
	for (int i = 0; i < numOfDevices; i++){
		safe_cuda(cudaSetDevice(i));
		
		for (int j = 0; j < 2; j++){
		
		dev_fields[i][j].setDevice(i);
		dev_mx[i][j].setDevice(i);
		dev_my[i][j].setDevice(i);
		dev_mz[i][j].setDevice(i);
		dev_total_mx[i][j].setDevice(i);
		dev_total_my[i][j].setDevice(i);
		dev_total_mz[i][j].setDevice(i);
		devStates[i][j].setDevice(i);
		dev_pars[i][j].setDevice(i);
		dev_basis[i][j].setDevice(i);
		
		dev_fields[i][j].malloc(meas_max);
		dev_mx[i][j].malloc(meas_max*number_of_particles);
		dev_my[i][j].malloc(meas_max*number_of_particles);
		dev_mz[i][j].malloc(meas_max*number_of_particles);
		dev_total_mx[i][j].malloc(meas_max);
		dev_total_my[i][j].malloc(meas_max);
		dev_total_mz[i][j].malloc(meas_max);
		devStates[i][j].malloc(number_of_particles);
		dev_pars[i][j].malloc();
		dev_basis[i][j] = basis;
		dev_basis[i][j].copyToDevice();
		cudaStreamCreate( &streams[i][j] );
		}
	}
	
	for (int devNum = 0; devNum < numOfDevices; devNum++){
	
	safe_cuda(cudaSetDevice(devNum));	
	int numOfSM = numOfSMPerDevice[devNum];
	
	for (int i = plan[devNum]; i < plan[devNum+1]; i=i+2){
	
		int steps = acqSet[i].getNumOfSteps();
		int steps1 = acqSet[i+1].getNumOfSteps();
		int measurements = acqSet[i].getNumOfMeas();
		int measurements1 = acqSet[i+1].getNumOfMeas();
		host_fields[i] = acqSet[i].getFieldFunctors();
		host_fields[i+1] = acqSet[i+1].getFieldFunctors();
		host_pars[i] = SimuParams(number_of_particles,
							   steps,
							   measurements, 
							   timestep, 
							   acqSet[i].getSeed(),
							   initialM.x, 
							   initialM.y, 
							   initialM.z ); 	


		host_pars[i+1] = SimuParams(number_of_particles,
							   steps1,
							   measurements1, 
							   timestep, 
							   acqSet[i+1].getSeed(),
							   initialM.x, 
							   initialM.y, 
							   initialM.z ); 	
		


	host_fields[i].copyToDevice(dev_fields[devNum][0], streams[devNum][0]);
	host_fields[i+1].copyToDevice(dev_fields[devNum][1], streams[devNum][1]);
		  
	  
	host_pars[i].copyToDevice(dev_pars[devNum][0], streams[devNum][0]);
	host_pars[i+1].copyToDevice(dev_pars[devNum][1], streams[devNum][1]);
	setup_kernel <<<blocks, threads, 0, streams[devNum][0]>>> (devStates[devNum][0].getPointer(), dev_pars[devNum][0].getPointer());
	setup_kernel <<<blocks, threads, 0, streams[devNum][1]>>> (devStates[devNum][1].getPointer(), dev_pars[devNum][1].getPointer());

	updateWalkersMag<<<blocks, threads, 0, streams[devNum][0]>>> (dev_pars[devNum][0].getPointer(),dev_basis[devNum][0].getPointer(), dev_fields[devNum][0].getPointer(), devStates[devNum][0].getPointer(), dev_mx[devNum][0].getPointer(), dev_my[devNum][0].getPointer(), dev_mz[devNum][0].getPointer() ) ;
	updateWalkersMag<<<blocks, threads, 0, streams[devNum][1]>>> (dev_pars[devNum][1].getPointer(),dev_basis[devNum][1].getPointer(), dev_fields[devNum][1].getPointer(), devStates[devNum][1].getPointer(), dev_mx[devNum][1].getPointer(), dev_my[devNum][1].getPointer(), dev_mz[devNum][1].getPointer() );

		dev_mz[devNum][0].sum(dev_total_mz[devNum][0], 768, numOfSM, measurements, number_of_particles, streams[devNum][0]);
		dev_mz[devNum][1].sum(dev_total_mz[devNum][1], 768, numOfSM, measurements1, number_of_particles, streams[devNum][1]);	
		dev_my[devNum][0].sum(dev_total_my[devNum][0], 768, numOfSM, measurements, number_of_particles, streams[devNum][0]);
		dev_my[devNum][1].sum(dev_total_my[devNum][1], 768, numOfSM, measurements1, number_of_particles, streams[devNum][1]);
		dev_mx[devNum][0].sum(dev_total_mx[devNum][0], 768, numOfSM, measurements, number_of_particles, streams[devNum][0]);
		dev_mx[devNum][1].sum(dev_total_mx[devNum][1], 768, numOfSM, measurements1, number_of_particles, streams[devNum][1]);
		
		host_mz[i].copyFromDevice(dev_total_mz[devNum][0], streams[devNum][0]);
		host_mz[i+1].copyFromDevice(dev_total_mz[devNum][1], streams[devNum][1]);	
		host_my[i].copyFromDevice(dev_total_my[devNum][0], streams[devNum][0]);
		host_my[i+1].copyFromDevice(dev_total_my[devNum][1], streams[devNum][1]);
		host_mx[i].copyFromDevice(dev_total_mx[devNum][0], streams[devNum][0]);
		host_mx[i+1].copyFromDevice(dev_total_mx[devNum][1], streams[devNum][1]);


				
    
		}
	
	}
	
	for (int devNum = 0; devNum < numOfDevices; devNum++){
		safe_cuda(cudaSetDevice(devNum));	
		cudaStreamSynchronize( streams[devNum][0] );
		cudaStreamSynchronize( streams[devNum][1] );
		cudaStreamDestroy( streams[devNum][0] );
		cudaStreamDestroy( streams[devNum][1] );

	}

	// ctimer.stop();
	// ctimer.display();
	
    addMagnetizationSet(host_mx, host_my, host_mz);
			
  }  
  
   void addMagnetizationSet(vector<pinnedVector<real> > & host_mx, vector<pinnedVector<real> > & host_my, vector<pinnedVector<real> > & host_mz){
    for (int i = 0; i < acqSet.size(); i++){
		host_mx[i].copyTo(acqSet[i].getMx());
		host_my[i].copyTo(acqSet[i].getMy());
		host_mz[i].copyTo(acqSet[i].getMz());
    }
  }
		
  void changeSeeds(long int seed){
	for(int i = 0; i < acqSet.size(); i++){
		acqSet[i].getSeed() = seed*(i+1);
	}
  }
  
  void calcEveryADC(){
      for (int i = 0; i < acqSet.size(); i++){
		getAcquisition(i).calcADC();
     } 	
  }
  
  
};
