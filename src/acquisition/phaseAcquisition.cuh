template <class FieldType> class phaseAcquisition {

  private:
    int numOfMeasurements;
    int numOfSteps;
    int lastAllocMeasurement;
    int seed;

    std::vector<real> mx;
    std::vector<real> my;
    std::vector<FieldType> gradientFunctors;
    std::vector<real> T2factor;
    real ADC;
    real ADCmeansquare;
    int points;

  public:
    phaseAcquisition(int _numOfMeas, int _numOfSteps, int _numOfParticles,
                     int _seed) {
        ADCmeansquare = 0;
        points = 0;
        mx = std::vector<real>(_numOfMeas, 0.0);
        my = std::vector<real>(_numOfMeas, 0.0);
        T2factor = std::vector<real>(_numOfParticles, 1.0);
        gradientFunctors = std::vector<FieldType>(_numOfMeas);
        numOfMeasurements = _numOfMeas;
        seed = _seed;
        numOfSteps = _numOfSteps;
        lastAllocMeasurement = 0;
    }

    phaseAcquisition() {}

    void addMeasurement(FieldType &field) {

        gradientFunctors[lastAllocMeasurement] = field;
        lastAllocMeasurement++;
    }

    std::vector<FieldType> &getGradientFunctors() { return gradientFunctors; }

    std::vector<real> &getMx() { return mx; }

    std::vector<real> &getMy() { return my; }

    real getSignal(int i) { return sqrt(mx[i] * mx[i] + my[i] * my[i]); }

    int getNumOfMeas() { return numOfMeasurements; }

    int getNumOfSteps() { return numOfSteps; }

    int &getSeed() { return seed; }

    void calcADC() {

        points++;
        std::vector<real> lnsignal;
        std::vector<real> bvalues;

        for (int i = 0; i < numOfMeasurements; i++) {

            bvalues.push_back(gradientFunctors[i].getB());
            lnsignal.push_back(log(sqrt(mx[i] * mx[i] + my[i] * my[i])));
        }

        if (points == 1) {
            ADC = -1.0 * linear_regression(lnsignal, bvalues);
            ADCmeansquare = ADC * ADC;
        }

        else {
            real ADCtemp = -1.0 * linear_regression(lnsignal, bvalues);
            ADC = (ADC * (points - 1) + ADCtemp) / (points);
            ADCmeansquare =
                (ADCmeansquare * (points - 1) + ADCtemp * ADCtemp) / (points);
        }
    }

    real &getADC() { return ADC; }

    real &getADCMeanSquare() { return ADCmeansquare; }

    int &getPoints() { return points; }

    real getADCError() { return sqrt(ADCmeansquare - ADC * ADC); }
};
