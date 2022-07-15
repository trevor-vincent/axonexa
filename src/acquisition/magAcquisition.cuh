template <class FieldType> class magAcquisition {

  private:
    int numOfMeasurements;
    int numOfSteps;
    int lastAllocMeasurement;
    int seed;

    std::vector<real> mx;
    std::vector<real> my;
    std::vector<real> mz;
    std::vector<FieldType> fieldFunctors;

    real ADC;
    real ADCmeansquare;
    int points;

  public:
    magAcquisition(int _numOfMeas, int _numOfSteps, int _seed) {
        points = 0;
        mx = std::vector<real>(_numOfMeas, 0.0);
        my = std::vector<real>(_numOfMeas, 0.0);
        mz = std::vector<real>(_numOfMeas, 0.0);
        fieldFunctors = std::vector<FieldType>(_numOfMeas);
        numOfMeasurements = _numOfMeas;
        seed = _seed;
        numOfSteps = _numOfSteps;
        lastAllocMeasurement = 0;
    }

    void addMeasurement(FieldType &field) {

        fieldFunctors[lastAllocMeasurement] = field;
        lastAllocMeasurement++;
    }

    std::vector<FieldType> &getFieldFunctors() { return fieldFunctors; }

    std::vector<real> &getMx() { return mx; }

    std::vector<real> &getMy() { return my; }

    std::vector<real> &getMz() { return mz; }

    int getNumOfMeas() { return numOfMeasurements; }

    int getNumOfSteps() { return numOfSteps; }

    int &getSeed() { return seed; }

    void calcADC() {

        points++;
        std::vector<real> lnsignal;
        std::vector<real> bvalues;

        for (int i = 0; i < numOfMeasurements; i++) {

            bvalues.push_back(fieldFunctors[i].getB());
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

    real getADC() { return ADC; }

    real getADCError() { return sqrt(ADCmeansquare - ADC * ADC); }
};
