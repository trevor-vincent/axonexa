// not finished

template <class T> class MeasurementDB {

  private:
    vector<T> _db;

  public:
    // Matrix<double> ADCtensor;
    double ADC;

    MeasurementDB() {}
    ~MeasurementDB() {}

    void add(T measurement) { _db.push_back(measurement); }

    void updateAll() {

        for (int i = 0; i < _db.size(); i++) {

            _db[i].updatePhase();
        }
    }

    void calcD() {

        vector<double> lnsignal;
        vector<double> b;
        for (int i = 0; i < _db.size(); i++) {
        }
    }
}