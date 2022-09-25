/*

This header file contains some very useful subroutines used throughout the code
which includes:

linear regression subroutines
multivariate linear regression subroutines
mean
std dev
mean ADC
FA
etc

*/

/* header files for the Numerical Recipes SVD fit */
#include "fitsvd.h"
#include "svd.h"
//#include "mrisequence.h"

/***************************************************************************************
        What: mean subrouutine
        Description:

        calculates the mean
 ****************************************************************************************/
double mean(double nums[], int number) {

    double mean = 0.0;
    for (int i = 0; i < number; i++) {
        mean = mean + nums[i];
    }
    mean = mean / number;
    return mean;
}

/***************************************************************************************
        What: mean_square subrouutine
        Description:

        calculates the mean of the square
 ****************************************************************************************/
double mean_square(double nums[], int number) {
    double mean = 0.0;
    for (int i = 0; i < number; i++) {
        mean = mean + nums[i] * nums[i];
    }
    mean = mean / number;
    return mean;
}

/***************************************************************************************
        What: stddev subrouutine
        Description:

        calculates the standard deviation
 ****************************************************************************************/
double std_dev(double nums[], int number) {

    double mean_num = mean(nums, number);
    double std_dev = 0.0;
    for (int i = 0; i < number; i++) {
        std_dev =
            std_dev + (nums[i] - mean_num) * (nums[i] - mean_num) / number;
    }
    std_dev = sqrt(std_dev);
    return std_dev;
}

/***************************************************************************************
        What: linear regression subroutine #1
        Description:

        Linear regression a la Bevington. He includes the uncertainty in yi
 (sigma_i) into the calculation of the slope. term 1 -> sum(1/sigma_i^2) term 2
 -> sum(yi*xi/sigma_i^2) term 3 -> sum(xi/sigma_i^2) term 4 -> sum(yi/sigma_i^2)
        term 5 -> sum(xi^2/sigma_i^2)

        returns slope
 ****************************************************************************************/

double linear_regression(double y[], double x[], int number, double sigma[],
                         double *error) {

    double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0, term5 = 0.0,
           delta;

    for (int i = 0; i < number; i++) {
        term1 = term1 + (1 / (sigma[i] * sigma[i]));
        term2 = term2 + (x[i] * y[i] / (sigma[i] * sigma[i]));
        term3 = term3 + (x[i] / (sigma[i] * sigma[i]));
        term4 = term4 + (y[i] / (sigma[i] * sigma[i]));
        term5 = term5 + (x[i] * x[i] / (sigma[i] * sigma[i]));
    }

    delta = term1 * term5 - term3 * term3;
    *error = sqrt((1 / delta) * term1);

    return (1 / delta) * (term1 * term2 - term3 * term4);
}

/***************************************************************************************
        What: linear regression subroutine #2
        Description:

        Linear regression without the uncertainty in y. This makes the vital
 assumption that the uncertainty in the yi's is CONSTANT. They should, however,
 be constant for us.
 ****************************************************************************************/

double linear_regression(double y[], double x[], int number) {

    double term1 = number, term2 = 0.0, term3 = 0.0, term4 = 0.0, term5 = 0.0,
           delta;

    for (int i = 0; i < number; i++) {
        term2 = term2 + (x[i] * y[i]);
        term3 = term3 + (x[i]);
        term4 = term4 + (y[i]);
        term5 = term5 + (x[i] * x[i]);
    }

    delta = term1 * term5 - term3 * term3;
    return (1 / delta) * (term1 * term2 - term3 * term4);
}

/***************************************************************************************
        What: linear regression subroutine #2
        Description:

        Linear regression without the uncertainty in y. This makes the vital
 assumption that the uncertainty in the yi's is CONSTANT. They should, however,
 be constant for us.
 ****************************************************************************************/

double linear_regression(vector<double> y, vector<double> x) {

    if (y.size() != x.size()) {
        std::cout << "ERROR LINEAR_REGRESSION Y.SIZE != X.SIZE " << std::endl;
        exit(0);
    }
    double term1 = y.size(), term2 = 0.0, term3 = 0.0, term4 = 0.0, term5 = 0.0,
           delta;

    for (int i = 0; i < y.size(); i++) {
        term2 = term2 + (x[i] * y[i]);
        term3 = term3 + (x[i]);
        term4 = term4 + (y[i]);
        term5 = term5 + (x[i] * x[i]);
    }

    delta = term1 * term5 - term3 * term3;
    return (1 / delta) * (term1 * term2 - term3 * term4);
}

/***************************************************************************************
        What: Fractional Anisotropy calculation subroutine
        Description:

        calculates FA
 ****************************************************************************************/
double FA_calc(double ADCx, double ADCy, double ADCz) {

    double mean_ADC, FA;
    mean_ADC = (ADCx + ADCy + ADCz) / 3;
    FA = sqrt(3.0 / 2.0) *
         sqrt((ADCx - mean_ADC) * (ADCx - mean_ADC) +
              (ADCy - mean_ADC) * (ADCy - mean_ADC) +
              (ADCz - mean_ADC) * (ADCz - mean_ADC)) /
         sqrt(ADCx * ADCx + ADCy * ADCy + ADCz * ADCz);
    return FA;
}

/***************************************************************************************
        What: multivariate linear regression subroutine #1
        Description:

        This defines the type of linear regression we want

        y(x1,x2,x3,x4,x5,x6) = a1 + a2*x1 + a3*x2 + a4*x3 + a5*x4 + a6*x5 +
 a7*x6

        the numerical recipes SVDFIT requires a function which defines y without
 the parameters a1->a7 and this is it.
 ****************************************************************************************/
VecDoub lin_multi(VecDoub_I &xx) {

    VecDoub ans(7);
    Doub x1 = xx[0];
    Doub x2 = xx[1];
    Doub x3 = xx[2];
    Doub x4 = xx[3];
    Doub x5 = xx[4];
    Doub x6 = xx[5];

    ans[0] = 1;
    ans[1] = x1;
    ans[2] = x2;
    ans[3] = x3;
    ans[4] = x4;
    ans[5] = x5;
    ans[6] = x6;

    return ans;
}

/***************************************************************************************
        What: multivariate linear regression subroutine #2
        Description:

        This is meat of the multivariate linear regression routine

 ****************************************************************************************/
double multi_lin(double lnsignal[], double b[][6], int num, double D[]) {

    MatDoub xx(num, 6); // n x 6 matrix where n is the number of b vectors and 6
                        // is the b components
    VecDoub yy(num);    // nx1 vector for the y values
    VecDoub sig(num);   // error

    for (int i = 0; i < num; i++) {

        // get y points
        yy[i] = lnsignal[i];
        sig[i] = 1.0;

        for (int j = 0; j < 6; j++) {

            // x points are the b vector points
            xx[i][j] = b[i][j];
        }
    }

    // runs the fit
    Fitsvd sol(xx, yy, sig, lin_multi, 1E-10);
    sol.fit();

    std::cout << "Tolerance = " << sol.tol << std::endl;
    std::cout << "Chi Squared = " << sol.chisq << std::endl;

    for (int i = 1; i < 7; i++) {

        // parameters the fit solves for (the D's)
        D[i - 1] = sol.a[i];
    }
}

/***************************************************************************************
        What: meanADC
        Description:

        Calculates the average and the std deviation of the adcs for back to
 back simulations

 ****************************************************************************************/
void meanADC(double ADC[], int num_of_inc, int num_of_repeat, double ADC_avg[],
             double ADC_error[]) {

    double *trial_ADC = new double[num_of_repeat];

    for (int i = 0; i < num_of_inc; i++) {

        for (int j = 0; j < num_of_repeat; j++) {

            trial_ADC[j] = ADC[j * num_of_inc + i];
        }

        ADC_avg[i] = mean(trial_ADC, num_of_repeat);
        ADC_error[i] = std_dev(trial_ADC, num_of_repeat);
    }

    delete[] trial_ADC;
}

void meanADC(const vector<double> &ADC, int num_of_inc, int num_of_repeat,
             vector<double> &ADC_avg, vector<double> &ADC_error) {

    double *trial_ADC = new double[num_of_repeat];

    for (int i = 0; i < num_of_inc; i++) {

        for (int j = 0; j < num_of_repeat; j++) {

            trial_ADC[j] = ADC[j * num_of_inc + i];
        }

        ADC_avg[i] = mean(trial_ADC, num_of_repeat);
        ADC_error[i] = std_dev(trial_ADC, num_of_repeat);
    }

    delete[] trial_ADC;
}

// template <class SEQUENCE>
// void meanADC (vector <MeasurementDB<SEQUENCE> > MDB_array, double & mean,
// double & stddev) {

// vector <double> ADC(MDB_array.size());
// mean = 0.0;
// double meansquare = 0.0;

// for ( int i = 0; i < MDB_array.size(); i++ ) {
// ADC[i] = MDB_array[i].calc_ADC();
// mean += ADC[i];
// meansquare += (ADC[i])*(ADC[i]);
// }

// mean /= MDB_array.size();
// meansquare /= MDB_array.size();

// stddev = sqrt(meansquare - mean);
// }

bool memberofarray(double num, double array[], int size_of_array) {

    for (int i = 0; i < size_of_array; i++) {

        if (doub_equal(num, array[i]) == true) {
            return true;
        }
    }

    return false;
}
