#include <iostream>
#include <vector>
#include <string>
#include <random>

using namespace std;


vector<double> arma_filter(vector<double> thetas, vector<double> phis, vector<double> x) {
    int phis_len = phis.size();
    int thetas_len = thetas.size();
    int x_len = x.size();
    vector<double> y(x_len, 0.0);

    for (int i = 0; i < x_len; i++) {
        double tmp = 0;
        // Moving Average
        for (int j = 1; j < thetas_len; j++) {
            if (i-j < 0) {
                continue;
            }
            tmp += thetas[j] * x[i-j];
        }

        // Autoregressive
        for (int k = 1; k < phis_len; k++) {
            if (i-k < 0) {
                continue;
            }
            tmp -= phis[k] * y[i-k];
        }
        y[i] = tmp;
    }
    return y;

}

// thetas: MA coefficients
// phis: AR coefficients
vector<double> arma2ma(vector<double> thetas, vector<double> phis, int order) {
    int t = thetas.size();
    int p = phis.size();
    vector<double> result(order, 0.0);
    vector<double> n_phis;

    // negate elements in phis
    for (int i = 0; i < p; i++) {
        n_phis.push_back(-phis[i]);
    }
    
    // padding
    if (t < order) {
        for (int i = t; i < order; i++) {
            thetas.push_back(0.0);
        }
    }
    if (p < order) {
        for (int i = p; i < order; i++) {
            n_phis.push_back(0.0);
        }
    }
    
    result[0] = 1;

    for (int i = 0; i < order-1; i++) {
        result[i+1] = thetas[i+1];

        for (int j = 0; j < (min(i, p)+1); j++) {
            // issue: multiplication levels of result[i+1] and result[i-j]*n_phis[j+1] must be the same
            // to allow addition
            // calculating result is a recursive process
            result[i+1] = result[i+1] + result[i-j]*n_phis[j+1];
        }
    }
    return result;
}

int main() {

    vector<double> maparams{0.35};
    vector<double> arparams{1.0, -0.65};
    int lags = 10;

    vector<double> arma2ma_coeff;

    arma2ma_coeff = arma2ma(maparams, arparams, lags);

    for (int i = 0; i < arma2ma_coeff.size(); i++) {
        cout << arma2ma_coeff[i] << endl;
    }

}