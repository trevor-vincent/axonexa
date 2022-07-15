// add devstates to kernel stream
template <class FieldType> class phaseAcquisitionStream {

  private:
    int steps_max;
    int meas_max;

    // put this into simuparams and remove it from here
    int number_of_particles;

  public:
    std::vector<phaseAcquisition<FieldType>> acqSet;

    phaseAcquisitionStream(int num) {
        number_of_particles = num;
        steps_max = -1;
        meas_max = -1;
    }
    phaseAcquisitionStream() {}
    ~phaseAcquisitionStream() {}

    void addAcquisition(phaseAcquisition<FieldType> acq) {
        acqSet.push_back(acq);
        if (acq.getNumOfSteps() > steps_max) {
            steps_max = acq.getNumOfSteps();
        }
        if (acq.getNumOfMeas() > meas_max) {
            meas_max = acq.getNumOfMeas();
        }
    }

    phaseAcquisition<FieldType> &getAcquisition(int num) { return acqSet[num]; }

    // put in Acquisition class
    template <class Basis>
    void runAcquisition(int num, Basis &basis, real timestep, int blocks,
                        int threads, int numOfSM) {

        int steps = acqSet[num].getNumOfSteps();
        int measurements = acqSet[num].getNumOfMeas();

        std::cout << "measurements = " << measurements << std::endl;
        std::cout << "number_of_particles = " << number_of_particles
                  << std::endl;

        cudaVector<FieldType> dev_fields(
            acqSet[num].getGradientFunctors().size());

        cudaVector<real> dev_phase(measurements * number_of_particles);

        cudaVector<curandState> dev_states(number_of_particles);
        cudaScalar<Basis> dev_basis(basis);
        cudaScalar<SimuParamsPhase> dev_par(
            SimuParamsPhase(number_of_particles, steps, measurements, timestep,
                            acqSet[num].getSeed(), (real)0.0));

        cudaVector<real> dev_total_mx(measurements);
        cudaVector<real> dev_total_my(measurements);
        cudaVector<real> dev_total_mz(measurements);

        dev_fields = acqSet[num].getGradientFunctors();
        dev_fields.copyToDevice();
        dev_basis.copyToDevice();
        dev_par.copyToDevice();
        setup_kernel<<<blocks, threads>>>(dev_states.getPointer(),
                                          dev_par.getPointer());

        updateWalkersPhase<FieldType, Basis><<<blocks, threads>>>(

            dev_par.getPointer(), dev_basis.getPointer(),
            dev_fields.getPointer(), dev_states.getPointer(),
            dev_phase.getPointer()

        );

        dev_phase.transformAndSum<COS_UNARY<real>>(
            dev_total_mx, 1024, numOfSM, measurements, number_of_particles, 0);
        dev_phase.transformAndSum<SIN_UNARY<real>>(
            dev_total_my, 1024, numOfSM, measurements, number_of_particles, 0);

        dev_total_mx.copyFromDevice();
        dev_total_my.copyFromDevice();

        dev_total_mx.copyTo(acqSet[num].getMx());
        dev_total_my.copyTo(acqSet[num].getMy());
    }

    template <class Basis>
    void runAcquisitionCPU(int num, Basis &basis, real timestep) {

        int steps = acqSet[num].getNumOfSteps();
        int measurements = acqSet[num].getNumOfMeas();

        std::vector<real> phase(number_of_particles * measurements);
        std::vector<real> total_mx(measurements);
        std::vector<real> total_my(measurements);

        SimuParamsPhase pars(number_of_particles, steps, measurements, timestep,
                             acqSet[num].getSeed(), (real)0.0);

        srand(acqSet[num].getSeed());

        updateWalkersPhaseCPU<FieldType, Basis>(
            &pars, &basis, &acqSet[num].getGradientFunctors()[0], &phase[0]);

        for (int i = 0; i < measurements; i++) {

            real total_mx_temp = 0.0;
            real total_my_temp = 0.0;

            for (int j = 0; j < number_of_particles; j++) {
                total_mx_temp += cos(phase[j + i * number_of_particles]);
                total_my_temp += sin(phase[j + i * number_of_particles]);
            }

            total_mx[i] = total_mx_temp;
            total_my[i] = total_my_temp;
        }

        acqSet[num].getMx() = total_mx;
        acqSet[num].getMy() = total_my;
    }

    template <class Basis>
    void runAcquisitionLattice(int num, Basis *basis, Lattice &lat,
                               real timestep, int blocks, int threads,
                               int numOfSM) {

        int steps = acqSet[num].getNumOfSteps();
        int measurements = acqSet[num].getNumOfMeas();

        cudaVector<FieldType> dev_fields(
            acqSet[num].getGradientFunctors().size());
        cudaVector<real> dev_phase(measurements * number_of_particles);
        cudaVector<real> dev_T2factors(number_of_particles);
        cudaVector<curandState> dev_states(number_of_particles);
        cudaVector<Basis> dev_basis(lat.getBasisSize());
        cudaScalar<Lattice> dev_lat(lat);

        for (int i = 0; i < dev_basis.size(); i++) {
            dev_basis[i] = basis[i];
        }

        cudaScalar<SimuParamsPhase> dev_par(
            SimuParamsPhase(number_of_particles, steps, measurements, timestep,
                            acqSet[num].getSeed(), (real)0.0));

        cudaVector<real> dev_total_mx(measurements);
        cudaVector<real> dev_total_my(measurements);
        cudaVector<real> dev_total_mz(measurements);

        dev_fields = acqSet[num].getGradientFunctors();
        dev_fields.copyToDevice();
        dev_basis.copyToDevice();
        dev_par.copyToDevice();
        dev_lat.copyToDevice();

        setup_kernel<<<blocks, threads>>>(dev_states.getPointer(),
                                          dev_par.getPointer());

        updateWalkersLatticePhase<FieldType, Basis><<<blocks, threads>>>(

            dev_par.getPointer(), dev_lat.getPointer(), dev_basis.getPointer(),
            dev_fields.getPointer(), dev_states.getPointer(),
            dev_T2factors.getPointer(), dev_phase.getPointer()

        );

#if defined USE_RELAXATION
        dev_phase.transformAndSumTwoVectors<COS_EXP_BINARY<real>>(
            dev_total_mx, dev_T2factors, 1024, numOfSM, measurements,
            number_of_particles, 0);
        dev_phase.transformAndSumTwoVectors<SIN_EXP_BINARY<real>>(
            dev_total_my, dev_T2factors, 1024, numOfSM, measurements,
            number_of_particles, 0);
#else
        dev_phase.transformAndSum<COS_UNARY<real>>(
            dev_total_mx, 1024, numOfSM, measurements, number_of_particles, 0);
        dev_phase.transformAndSum<SIN_UNARY<real>>(
            dev_total_my, 1024, numOfSM, measurements, number_of_particles, 0);
#endif

        dev_total_mx.copyFromDevice();
        dev_total_my.copyFromDevice();

        dev_total_mx.copyTo(acqSet[num].getMx());
        dev_total_my.copyTo(acqSet[num].getMy());
    }

    template <class Basis>
    void runAcquisitionWC(int num, Basis &basis, real timestep, int blocks,
                          int threads) {

        int steps = acqSet[num].getNumOfSteps();
        int measurements = acqSet[num].getNumOfMeas();

        cudaVector<real> dev_x(number_of_particles);
        cudaVector<real> dev_y(number_of_particles);
        cudaVector<real> dev_z(number_of_particles);
        cudaVector<real> dev_Gx(steps);
        cudaVector<real> dev_Gy(steps);
        cudaVector<real> dev_Gz(steps);

        cudaVector<real> dev_phase(number_of_particles);
        cudaVector<unsigned int> devStates(number_of_particles);
        cudaScalar<Basis> dev_basis(basis);
        cudaScalar<SimuParamsWC> dev_par(SimuParamsWC(steps, timestep));
        dev_basis.copyToDevice();
        dev_par.copyToDevice();

        for (int i = 0; i < measurements; i++) {

            for (int j = 0; j < steps; j++) {
                dev_Gx[j] =
                    acqSet[num].getGradientFunctors()[i](j * timestep).x;
                dev_Gy[j] =
                    acqSet[num].getGradientFunctors()[i](j * timestep).y;
                dev_Gz[j] =
                    acqSet[num].getGradientFunctors()[i](j * timestep).z;
            }

            srand(acqSet[num].getSeed());
            for (int j = 0; j < number_of_particles; j++) {
                devStates[j] = (unsigned int)rand();
                Vector3 randVec = basis.unifRandCPU();
                dev_x[j] = randVec.x;
                dev_y[j] = randVec.y;
                dev_z[j] = randVec.z;
            }

            dev_x.copyToDevice();
            dev_y.copyToDevice();
            dev_z.copyToDevice();
            dev_Gx.copyToDevice();
            dev_Gy.copyToDevice();
            dev_Gz.copyToDevice();
            devStates.copyToDevice();

            kernelWC<<<blocks, threads>>>(
                devStates.getPointer(), dev_basis.getPointer(),
                dev_par.getPointer(), dev_x.getPointer(), dev_y.getPointer(),
                dev_z.getPointer(), dev_phase.getPointer(), dev_Gx.getPointer(),
                dev_Gy.getPointer(), dev_Gz.getPointer());

            acqSet[num].getMx()[i] =
                computePhaseSx(number_of_particles, dev_phase.getPointer(), 1);
            acqSet[num].getMy()[i] =
                computePhaseSx(number_of_particles, dev_phase.getPointer(), 0);
        }
    }

    template <class Basis>
    void runAcquisitionStream(Basis &basis, real timestep, int blocks,
                              int threads, int numOfDevices,
                              std::vector<int> &plan,
                              std::vector<int> &numOfSMPerDevice) {

        vector<pinnedVector<FieldType>> host_fields(acqSet.size());
        vector<pinnedScalar<SimuParamsPhase>> host_pars(acqSet.size());
        vector<pinnedVector<real>> host_mx(acqSet.size());
        vector<pinnedVector<real>> host_my(acqSet.size());

        for (int i = 0; i < acqSet.size(); i++) {

            host_fields[i].alloc(meas_max);
            host_mx[i].alloc(meas_max);
            host_my[i].alloc(meas_max);
        }

        std::vector<std::vector<cudaStream_t>> streams(
            numOfDevices, std::vector<cudaStream_t>(2));
        std::vector<std::vector<cudaVector<curandState>>> devStates(
            numOfDevices, std::vector<cudaVector<curandState>>(2));

        std::vector<std::vector<cudaVector<FieldType>>> dev_fields(
            numOfDevices, std::vector<cudaVector<FieldType>>(2));
        std::vector<std::vector<cudaScalar<SimuParamsPhase>>> dev_pars(
            numOfDevices, std::vector<cudaScalar<SimuParamsPhase>>(2));
        std::vector<std::vector<cudaScalar<Basis>>> dev_basis(
            numOfDevices, std::vector<cudaScalar<Basis>>(2));

        std::vector<std::vector<cudaVector<real>>> dev_phase(
            numOfDevices, std::vector<cudaVector<real>>(2));
        std::vector<std::vector<cudaVector<real>>> dev_total_mx(
            numOfDevices, std::vector<cudaVector<real>>(2));
        std::vector<std::vector<cudaVector<real>>> dev_total_my(
            numOfDevices, std::vector<cudaVector<real>>(2));

        for (int i = 0; i < numOfDevices; i++) {

            safe_cuda(cudaSetDevice(i));
            for (int j = 0; j < 2; j++) {

                dev_fields[i][j].setDevice(i);
                dev_phase[i][j].setDevice(i);
                dev_total_mx[i][j].setDevice(i);
                dev_total_my[i][j].setDevice(i);
                devStates[i][j].setDevice(i);
                dev_pars[i][j].setDevice(i);
                dev_basis[i][j].setDevice(i);

                dev_fields[i][j].malloc(meas_max);
                dev_phase[i][j].malloc(meas_max * number_of_particles);
                dev_total_mx[i][j].malloc(meas_max);
                dev_total_my[i][j].malloc(meas_max);
                devStates[i][j].malloc(number_of_particles);
                dev_pars[i][j].malloc();
                dev_basis[i][j] = basis;
                dev_basis[i][j].copyToDevice();
                cudaStreamCreate(&streams[i][j]);
            }
        }

        for (int devNum = 0; devNum < numOfDevices; devNum++) {

            safe_cuda(cudaSetDevice(devNum));
            int numOfSM = numOfSMPerDevice[devNum];

            for (int i = plan[devNum]; i < plan[devNum + 1]; i = i + 2) {

                int steps = acqSet[i].getNumOfSteps();
                int steps1 = acqSet[i + 1].getNumOfSteps();
                int measurements = acqSet[i].getNumOfMeas();
                int measurements1 = acqSet[i + 1].getNumOfMeas();
                host_fields[i] = acqSet[i].getGradientFunctors();
                host_fields[i + 1] = acqSet[i + 1].getGradientFunctors();
                host_pars[i] =
                    SimuParamsPhase(number_of_particles, steps, measurements,
                                    timestep, acqSet[i].getSeed(), (real)0.0);

                host_pars[i + 1] = SimuParamsPhase(
                    number_of_particles, steps1, measurements1, timestep,
                    acqSet[i + 1].getSeed(), (real)0.0);

                host_fields[i].copyToDevice(dev_fields[devNum][0],
                                            streams[devNum][0]);
                host_fields[i + 1].copyToDevice(dev_fields[devNum][1],
                                                streams[devNum][1]);
                host_pars[i].copyToDevice(dev_pars[devNum][0],
                                          streams[devNum][0]);
                host_pars[i + 1].copyToDevice(dev_pars[devNum][1],
                                              streams[devNum][1]);
                cudaStreamQuery(streams[devNum][0]);
                cudaStreamQuery(streams[devNum][1]);
                setup_kernel<<<blocks, threads, 0, streams[devNum][0]>>>(
                    devStates[devNum][0].getPointer(),
                    dev_pars[devNum][0].getPointer());
                setup_kernel<<<blocks, threads, 0, streams[devNum][1]>>>(
                    devStates[devNum][1].getPointer(),
                    dev_pars[devNum][1].getPointer());
                cudaStreamQuery(streams[devNum][0]);
                cudaStreamQuery(streams[devNum][1]);
                updateWalkersPhase<<<blocks, threads, 0, streams[devNum][0]>>>(
                    dev_pars[devNum][0].getPointer(),
                    dev_basis[devNum][0].getPointer(),
                    dev_fields[devNum][0].getPointer(),
                    devStates[devNum][0].getPointer(),
                    dev_phase[devNum][0].getPointer());
                updateWalkersPhase<<<blocks, threads, 0, streams[devNum][1]>>>(
                    dev_pars[devNum][1].getPointer(),
                    dev_basis[devNum][1].getPointer(),
                    dev_fields[devNum][1].getPointer(),
                    devStates[devNum][1].getPointer(),
                    dev_phase[devNum][1].getPointer());
                cudaStreamQuery(streams[devNum][0]);
                cudaStreamQuery(streams[devNum][1]);
                // #if defined USE_RELAXATION
                // dev_phase[devNum][0].transformAndSumTwoVectors<SIN_EXP_BINARY<real>
                // >(dev_total_my[devNum][0], dev_T2factors[devNum][0], 1024,
                // numOfSM, measurements, number_of_particles,
                // streams[devNum][0]);
                // dev_phase[devNum][1].transformAndSumTwoVectors<SIN_EXP_BINARY<real>
                // >(dev_total_my[devNum][1], dev_T2factors[devNum][1], 1024,
                // numOfSM, measurements, number_of_particles,
                // streams[devNum][1]);
                // dev_phase[devNum][0].transformAndSumTwoVectors<COS_EXP_BINARY<real>
                // >(dev_total_mx[devNum][0], dev_T2factors[devNum][0],  1024,
                // numOfSM, measurements, number_of_particles,
                // streams[devNum][0]);
                // dev_phase[devNum][1].transformAndSumTwoVectors<COS_EXP_BINARY<real>>(dev_total_mx[devNum][1],
                // dev_T2factors[devNum][1], 1024, numOfSM, measurements,
                // number_of_particles, streams[devNum][1]); #else
                dev_phase[devNum][0].transformAndSum<SIN_UNARY<real>>(
                    dev_total_my[devNum][0], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][0]);
                dev_phase[devNum][1].transformAndSum<SIN_UNARY<real>>(
                    dev_total_my[devNum][1], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][1]);
                dev_phase[devNum][0].transformAndSum<COS_UNARY<real>>(
                    dev_total_mx[devNum][0], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][0]);
                dev_phase[devNum][1].transformAndSum<COS_UNARY<real>>(
                    dev_total_mx[devNum][1], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][1]);
                // #endif
                cudaStreamQuery(streams[devNum][0]);
                cudaStreamQuery(streams[devNum][1]);

                host_my[i].copyFromDevice(dev_total_my[devNum][0],
                                          streams[devNum][0]);
                host_my[i + 1].copyFromDevice(dev_total_my[devNum][1],
                                              streams[devNum][1]);
                host_mx[i].copyFromDevice(dev_total_mx[devNum][0],
                                          streams[devNum][0]);
                host_mx[i + 1].copyFromDevice(dev_total_mx[devNum][1],
                                              streams[devNum][1]);
                cudaStreamQuery(streams[devNum][0]);
                cudaStreamQuery(streams[devNum][1]);
            }
        }

        for (int devNum = 0; devNum < numOfDevices; devNum++) {
            safe_cuda(cudaSetDevice(devNum));
            cudaStreamSynchronize(streams[devNum][0]);
            cudaStreamSynchronize(streams[devNum][1]);
            cudaStreamDestroy(streams[devNum][0]);
            cudaStreamDestroy(streams[devNum][1]);
        }

        addMagnetizationSet(host_mx, host_my);
    }

    template <class Basis>
    void runAcquisitionStreamLattice(Basis *basis, Lattice &lat, real timestep,
                                     int blocks, int threads, int numOfDevices,
                                     std::vector<int> &plan,
                                     std::vector<int> &numOfSMPerDevice) {

        vector<pinnedVector<FieldType>> host_fields(acqSet.size());
        vector<pinnedScalar<SimuParamsPhase>> host_pars(acqSet.size());
        vector<pinnedVector<real>> host_mx(acqSet.size());
        vector<pinnedVector<real>> host_my(acqSet.size());

        for (int i = 0; i < acqSet.size(); i++) {
            host_fields[i].alloc(meas_max);
            host_mx[i].alloc(meas_max);
            host_my[i].alloc(meas_max);
        }

        std::vector<std::vector<cudaStream_t>> streams(
            numOfDevices, std::vector<cudaStream_t>(2));
        std::vector<std::vector<cudaVector<curandState>>> devStates(
            numOfDevices, std::vector<cudaVector<curandState>>(2));

        std::vector<std::vector<cudaVector<FieldType>>> dev_fields(
            numOfDevices, std::vector<cudaVector<FieldType>>(2));
        std::vector<std::vector<cudaScalar<SimuParamsPhase>>> dev_pars(
            numOfDevices, std::vector<cudaScalar<SimuParamsPhase>>(2));
        std::vector<std::vector<cudaScalar<Lattice>>> dev_lat(
            numOfDevices, std::vector<cudaScalar<Lattice>>(2));
        std::vector<std::vector<cudaVector<Basis>>> dev_basis(
            numOfDevices, std::vector<cudaVector<Basis>>(2));

        std::vector<std::vector<cudaVector<real>>> dev_phase(
            numOfDevices, std::vector<cudaVector<real>>(2));
        std::vector<std::vector<cudaVector<real>>> dev_T2factors(
            numOfDevices, std::vector<cudaVector<real>>(2));
        std::vector<std::vector<cudaVector<real>>> dev_total_mx(
            numOfDevices, std::vector<cudaVector<real>>(2));
        std::vector<std::vector<cudaVector<real>>> dev_total_my(
            numOfDevices, std::vector<cudaVector<real>>(2));

        for (int i = 0; i < numOfDevices; i++) {

            safe_cuda(cudaSetDevice(i));
            for (int j = 0; j < 2; j++) {

                dev_fields[i][j].setDevice(i);
                dev_phase[i][j].setDevice(i);
                dev_T2factors[i][j].setDevice(i);
                dev_total_mx[i][j].setDevice(i);
                dev_total_my[i][j].setDevice(i);
                devStates[i][j].setDevice(i);
                dev_pars[i][j].setDevice(i);
                dev_lat[i][j].setDevice(i);
                dev_basis[i][j].setDevice(i);

                dev_fields[i][j].malloc(meas_max);
                dev_basis[i][j].malloc(lat.getBasisSize());
                dev_phase[i][j].malloc(meas_max * number_of_particles);
                dev_T2factors[i][j].malloc(number_of_particles);
                dev_total_mx[i][j].malloc(meas_max);
                dev_total_my[i][j].malloc(meas_max);
                devStates[i][j].malloc(number_of_particles);
                dev_pars[i][j].malloc();
                dev_lat[i][j].malloc();

                for (int k = 0; k < lat.getBasisSize(); k++) {
                    dev_basis[i][j][k] = basis[k];
                }

                dev_lat[i][j] = lat;
                dev_basis[i][j].copyToDevice();
                dev_lat[i][j].copyToDevice();
                cudaStreamCreate(&streams[i][j]);
            }
        }

        for (int devNum = 0; devNum < numOfDevices; devNum++) {

            safe_cuda(cudaSetDevice(devNum));
            int numOfSM = numOfSMPerDevice[devNum];

            for (int i = plan[devNum]; i < plan[devNum + 1]; i = i + 2) {

                int steps = acqSet[i].getNumOfSteps();
                int steps1 = acqSet[i + 1].getNumOfSteps();
                int measurements = acqSet[i].getNumOfMeas();
                int measurements1 = acqSet[i + 1].getNumOfMeas();
                host_fields[i] = acqSet[i].getGradientFunctors();
                host_fields[i + 1] = acqSet[i + 1].getGradientFunctors();
                host_pars[i] =
                    SimuParamsPhase(number_of_particles, steps, measurements,
                                    timestep, acqSet[i].getSeed(), (real)0.0);

                host_pars[i + 1] = SimuParamsPhase(
                    number_of_particles, steps1, measurements1, timestep,
                    acqSet[i + 1].getSeed(), (real)0.0);

                host_fields[i].copyToDevice(dev_fields[devNum][0],
                                            streams[devNum][0]);
                host_fields[i + 1].copyToDevice(dev_fields[devNum][1],
                                                streams[devNum][1]);

                host_pars[i].copyToDevice(dev_pars[devNum][0],
                                          streams[devNum][0]);
                host_pars[i + 1].copyToDevice(dev_pars[devNum][1],
                                              streams[devNum][1]);
                setup_kernel<<<blocks, threads, 0, streams[devNum][0]>>>(
                    devStates[devNum][0].getPointer(),
                    dev_pars[devNum][0].getPointer());
                setup_kernel<<<blocks, threads, 0, streams[devNum][1]>>>(
                    devStates[devNum][1].getPointer(),
                    dev_pars[devNum][1].getPointer());

                updateWalkersLatticePhase<<<blocks, threads, 0,
                                            streams[devNum][0]>>>(
                    dev_pars[devNum][0].getPointer(),
                    dev_lat[devNum][0].getPointer(),
                    dev_basis[devNum][0].getPointer(),
                    dev_fields[devNum][0].getPointer(),
                    devStates[devNum][0].getPointer(),
                    dev_T2factors[devNum][0].getPointer(),
                    dev_phase[devNum][0].getPointer());
                updateWalkersLatticePhase<<<blocks, threads, 0,
                                            streams[devNum][1]>>>(
                    dev_pars[devNum][1].getPointer(),
                    dev_lat[devNum][1].getPointer(),
                    dev_basis[devNum][1].getPointer(),
                    dev_fields[devNum][1].getPointer(),
                    devStates[devNum][1].getPointer(),
                    dev_T2factors[devNum][1].getPointer(),
                    dev_phase[devNum][1].getPointer());

#if defined USE_RELAXATION
                dev_phase[devNum][0]
                    .transformAndSumTwoVectors<SIN_EXP_BINARY<real>>(
                        dev_total_my[devNum][0], dev_T2factors[devNum][0], 1024,
                        numOfSM, measurements, number_of_particles,
                        streams[devNum][0]);
                dev_phase[devNum][1]
                    .transformAndSumTwoVectors<SIN_EXP_BINARY<real>>(
                        dev_total_my[devNum][1], dev_T2factors[devNum][1], 1024,
                        numOfSM, measurements, number_of_particles,
                        streams[devNum][1]);
                dev_phase[devNum][0]
                    .transformAndSumTwoVectors<COS_EXP_BINARY<real>>(
                        dev_total_mx[devNum][0], dev_T2factors[devNum][0], 1024,
                        numOfSM, measurements, number_of_particles,
                        streams[devNum][0]);
                dev_phase[devNum][1]
                    .transformAndSumTwoVectors<COS_EXP_BINARY<real>>(
                        dev_total_mx[devNum][1], dev_T2factors[devNum][1], 1024,
                        numOfSM, measurements, number_of_particles,
                        streams[devNum][1]);
#else
                dev_phase[devNum][0].transformAndSum<SIN_UNARY<real>>(
                    dev_total_my[devNum][0], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][0]);
                dev_phase[devNum][1].transformAndSum<SIN_UNARY<real>>(
                    dev_total_my[devNum][1], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][1]);
                dev_phase[devNum][0].transformAndSum<COS_UNARY<real>>(
                    dev_total_mx[devNum][0], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][0]);
                dev_phase[devNum][1].transformAndSum<COS_UNARY<real>>(
                    dev_total_mx[devNum][1], 1024, numOfSM, measurements,
                    number_of_particles, streams[devNum][1]);
#endif
                host_my[i].copyFromDevice(dev_total_my[devNum][0],
                                          streams[devNum][0]);
                host_my[i + 1].copyFromDevice(dev_total_my[devNum][1],
                                              streams[devNum][1]);
                host_mx[i].copyFromDevice(dev_total_mx[devNum][0],
                                          streams[devNum][0]);
                host_mx[i + 1].copyFromDevice(dev_total_mx[devNum][1],
                                              streams[devNum][1]);
            }
        }

        for (int devNum = 0; devNum < numOfDevices; devNum++) {
            safe_cuda(
                cudaSetDevice(devNum)); // ADDED THIS RECENTLY (MARCH 4th 2013)
            cudaStreamSynchronize(streams[devNum][0]);
            cudaStreamSynchronize(streams[devNum][1]);
            cudaStreamDestroy(streams[devNum][0]);
            cudaStreamDestroy(streams[devNum][1]);
        }

        addMagnetizationSet(host_mx, host_my);
    }

    void addMagnetizationSet(vector<pinnedVector<real>> &host_mx,
                             vector<pinnedVector<real>> &host_my) {
        for (int i = 0; i < acqSet.size(); i++) {
            host_mx[i].copyTo(acqSet[i].getMx());
            host_my[i].copyTo(acqSet[i].getMy());
        }
    }

    void changeSeeds(long int seed) {
        for (int i = 0; i < acqSet.size(); i++) {
            acqSet[i].getSeed() = seed * (i + 1);
        }
    }

    void calcEveryADC() {
        for (int i = 0; i < acqSet.size(); i++) {
            getAcquisition(i).calcADC();
        }
    }

    void flushADC() {
        for (int i = 0; i < acqSet.size(); i++) {
            getAcquisition(i).getADC() = 0.0;
            getAcquisition(i).getADCMeanSquare() = 0.0;
            getAcquisition(i).getPoints() = 0;
        }
    }
};
