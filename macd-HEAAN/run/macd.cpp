/*
* Moving Average Convergence Divergence with HEAAN Library
* The results are output into the data folder
* Run macd-heaan.ipynb to visualise the results
*/
#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <stdlib.h>
#include "../src/HEAAN.h"

using namespace std;

// ***************************************************************************************************
// Helper functions
// ***************************************************************************************************

// Import data in csv file into a vector of doubles
inline vector<double> csv2vec(string file) {
    vector<double> data;
    ifstream inputFile(file);
    long l = 0;

    // Iteratively read lines and push them into the vector
    while (inputFile) {
        l++;
        string line;
        if (!getline(inputFile, line)) {
            break;
        }
        try {
            data.push_back(stof(line));
        }
        catch (const std::invalid_argument e) {
            cout << "NaN found in file " << file << " line " << l
                 << endl;
            e.what();
        }
    }
    if (!inputFile.eof()) {
        cerr << "Could not read file " << file << "\n";
        __throw_invalid_argument("File not found.");
    }
    return data;
}

// Print a vector of doubles
inline void printVector(vector<double>& vec) {
    for(long i=0; i < vec.size(); i++) {
        cout << vec[i] << ' ';
    }
    cout << "\n";
}

// Return the sum of all the Ciphertexts in a vector
inline Ciphertext sum(vector<Ciphertext>& vec, Scheme &scheme) {
    Ciphertext result_encrypted = vec[0];
    for(long i = 1; i < vec.size(); i++) {
        Ciphertext val = vec[i];
        scheme.addAndEqual(result_encrypted, val);
    }
    return result_encrypted;
}

// Return a subset of the vector of Ciphertexts
inline vector<Ciphertext> slice(vector<Ciphertext>& vec, long start , long end) {
    long oldlen = vec.size();
    long newlen;
    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }
    vector<Ciphertext> res;
    for (long i = 0; i < newlen; i++) {
        res.push_back(vec[start+i]);
    }
    return res;
}

// ***************************************************************************************************
// MACD analysis functions
// ***************************************************************************************************

// Return the weighted moving averages of the input Ciphertexts in a moving window
inline vector<Ciphertext> wma(vector<Ciphertext>& data, long m, Scheme &scheme, long logq, long logp, long logn) {
    vector<Ciphertext> wma;
    vector<double> weights;

    // Generate a list of weights of the window size
    for (long i = 0; i < m; i++) {
        double w = (2.0*(i+1.0)/(m*(m+1.0)));
        weights.push_back(w);
    }

    // Multiply the corresponding Ciphertext and weight in a window and sum the results up 
    for (long i = 0; i < data.size()-m; i++) {
        vector<Ciphertext> data_sliced;
        vector<Ciphertext> window;

        // Get the data in the moving window
        data_sliced = slice(data, i, i+m);

        for (long j = 0; j < m; j++) {
            complex<double> tmp_weight;
			Ciphertext tmp = data_sliced[j];
			tmp_weight.real(weights[j]);

            // Multiply the Ciphertext and weight and rescale the result
            scheme.multByConstAndEqual(tmp, tmp_weight, logp);
			scheme.reScaleByAndEqual(tmp, logp);

            window.push_back(tmp);
        }

        // Sum the multiplication results up to get the weighted moving average
        Ciphertext res = sum(window, scheme);
        wma.push_back(res);
    }
    return wma;
}

// Return the trading decisions based on the MACD signals
// Denote MACD signals as m(t) where t is the time and sign function as sign()
// Decision(t) = (m(t-1)-m(t))*(sign(m(t-1)*m(t))-1)
// sign(m(t-1)*m(t))-1 indicates the turning point, where MACD signals cross the X axis
// m(t-1)-m(t) indicates the trend
// Decision(t) = 0 means "Do nothing" or "Hold"
// Decision(t) > 0 means "Buy"
// Decision(t) < 0 means "Sell"
// We use a polynomial approximation of tanh to approximate sign function
// The approximation of sign function: f(x) = 0.375*x^5 - 1.25*x^3 + 1.875*x where x = m(t)*m(t-1)
inline vector<Ciphertext> decision(vector<Ciphertext>& macd, Scheme &scheme, long logq, long logp) {
    vector<Ciphertext> decisions;
    for (int i = 1; i < macd.size(); i++) {

        // Store the latest two MACD signals m(t) and m(t-1)
        Ciphertext mt = macd[i];
        Ciphertext mt_1 = macd[i-1];

        // Set up coefficients
        complex<double> coeff1, coeff2, coeff3, offset, coeff_norm;
        coeff1.real(0.375);
        coeff2.real(-1.25);
        coeff3.real(1.875);
        offset.real(-1.0);
		coeff_norm.real(0.25);

        // Normalise the MACD signals so that most of the data is in range [-1, 1]
        scheme.multByConstAndEqual(mt, coeff_norm, logp);
        scheme.reScaleByAndEqual(mt, logp);
        scheme.multByConstAndEqual(mt_1, coeff_norm, logp);
        scheme.reScaleByAndEqual(mt_1, logp);

        // Calculate m(t-1)-m(t) and m(t-1)*m(t)
        // m(t-1)*m(t) is denoted as x in the following steps
        Ciphertext diff_res;
        scheme.sub(diff_res, mt_1, mt);
        Ciphertext mult_res;
        scheme.mult(mult_res, mt_1, mt);
        scheme.reScaleByAndEqual(mult_res, logp);

        // Calculate x^2
        Ciphertext x2_encrypted;
        scheme.square(x2_encrypted, mult_res);
        scheme.reScaleByAndEqual(x2_encrypted, logp);

        // Calculate x^3
        Ciphertext x3_encrypted;
        scheme.modDownToAndEqual(mult_res, x2_encrypted.logq);
        scheme.mult(x3_encrypted, mult_res, x2_encrypted);
        scheme.reScaleByAndEqual(x3_encrypted, logp);

        // Calculate x^5
        Ciphertext x5_encrypted;
        scheme.modDownToAndEqual(x2_encrypted, x3_encrypted.logq);
        scheme.mult(x5_encrypted, x2_encrypted, x3_encrypted);
        scheme.reScaleByAndEqual(x5_encrypted, logp);

        // Calculate 0.375*x^5
        scheme.multByConstAndEqual(x5_encrypted, coeff1, logp);
        scheme.reScaleByAndEqual(x5_encrypted, logp);

        // Calculate -1.25*x^3
        scheme.multByConstAndEqual(x3_encrypted, coeff2, logp);
        scheme.reScaleByAndEqual(x3_encrypted, logp);

        // Calculate 1.875*x
        scheme.multByConstAndEqual(mult_res, coeff3, logp);
        scheme.reScaleByAndEqual(mult_res, logp);

        // Calculate 0.375*x^5-1.25*x^3+1.875*x-1
        Ciphertext sign_encrypted;
        scheme.modDownToAndEqual(x3_encrypted, x5_encrypted.logq);
        scheme.add(sign_encrypted, x5_encrypted, x3_encrypted);
        scheme.modDownToAndEqual(mult_res, sign_encrypted.logq);
        scheme.addAndEqual(sign_encrypted, mult_res);
        scheme.addConstAndEqual(sign_encrypted, offset, logp);

        // Calculate (m(t-1)-m(t))*(f(x)-1)
        Ciphertext result_encrypted;
        scheme.modDownToAndEqual(diff_res, sign_encrypted.logq);
        scheme.mult(result_encrypted, sign_encrypted, diff_res);
        scheme.reScaleByAndEqual(result_encrypted, logp);

        decisions.push_back(result_encrypted);
    }
    return decisions;
}

// ***************************************************************************************************
// Data input/output functions
// ***************************************************************************************************

// Simulate the data import process. Return corresponding sample data according to the time
inline Ciphertext getSample(long time, long logq, long logp, long logn, Scheme& scheme) {
    // Read data from csv file
    vector<double> prices;
	prices = csv2vec("../data/apple_prices.csv");

    long n = pow(2, logn);
    complex<double>* sample = new complex<double>[n];
    if (time <= prices.size()) { 
        sample[0].real(prices[time]);
        for (long i = 1; i < n; ++i) {
            sample[i].real(0.0);
        }
    } else {
        cout << "Error: time exceeds sample data range! " << endl;
        for (long i = 0; i < n; ++i) {
            sample[i].real(0.0);
        }
    }
    Ciphertext sample_encrypted;
    scheme.encrypt(sample_encrypted, sample, n, logp, logq);
    return sample_encrypted;
}

// Append new sample data to the history data
inline void assembleSample(Ciphertext& sample, vector<Ciphertext>& past_prices) {
    past_prices.push_back(sample);
}

// Decrytp a vector of Ciphertexts into a vector of doubles
inline vector<double> decryptVec(vector<Ciphertext>& ct, Scheme& scheme, SecretKey& secretKey) {
    vector<double> res;
    for (int i = 0; i < ct.size(); i++) {
        complex<double>* val;
        val = scheme.decrypt(secretKey, ct[i]);
        res.push_back(real(val[0]));
    }
    return res;
}

// Write data to csv files for plotting
inline void exportData(vector<double> vec, string file) {
    ofstream output_file(file);
    ostream_iterator<double> output_iterator(output_file, "\n");
    copy(vec.begin(), vec.end(), output_iterator);
}

//********************************************************************************
// Main function
//********************************************************************************
int main() {

	long logq = 600; // Ciphertext polynomial modulus
	long logp = 30; // Scale
	long logn = 10; // n = number of slots in a Ciphertext
	long time_max = 200; // Sample data size

	SetNumThreads(8);
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	cout << "Data Import Begins" << endl;
	cout << "\n";

    vector<Ciphertext> prices_encrypted;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample(i, logq, logp, logn, scheme);
        assembleSample(sample, prices_encrypted);   
    }

	cout << "Data Import Finished" << endl;
	cout << "\n";

	cout << "MACD Analysis Begins" << endl;
	cout << "\n";

	vector<Ciphertext> wma12_encrypted = wma(prices_encrypted, 12, scheme, logq, logp, logn);
	vector<Ciphertext> wma12_encrypted_sliced = slice(wma12_encrypted, 14, wma12_encrypted.size());

	vector<Ciphertext> wma26_encrypted = wma(prices_encrypted, 26, scheme, logq, logp, logn);

	vector<Ciphertext> wma_diff_encrypted;
    for (int i = 0; i < wma26_encrypted.size(); i++) {
        Ciphertext tmp_diff;
        scheme.sub(tmp_diff, wma12_encrypted_sliced[i], wma26_encrypted[i]);
        wma_diff_encrypted.push_back(tmp_diff);
    }

    vector<Ciphertext> wma9_encrypted = wma(wma_diff_encrypted, 9, scheme, logq, logp, logn);

    vector<Ciphertext> wma_diff_sliced = slice(wma_diff_encrypted, 9, wma_diff_encrypted.size());

    vector<Ciphertext> macd_encrypted;
    for (int i = 0; i < wma9_encrypted.size(); i++) {
        Ciphertext tmp_macd;
        scheme.modDownToAndEqual(wma_diff_sliced[i], wma9_encrypted[i].logq);
        scheme.sub(tmp_macd, wma_diff_sliced[i], wma9_encrypted[i]);
        macd_encrypted.push_back(tmp_macd);
    }

	cout << "MACD Analysis Finished" << endl;
	cout << "\n";

    cout << "Decision Analysis Begins" << endl;
	cout << "\n";

	vector<Ciphertext> decisions_encrypted = decision(macd_encrypted, scheme, logq, logp);

    cout << "Decision Analysis Finished" << endl;
    cout << "\n";

    cout << "Data Export Begins" << endl;
	cout << "\n";
    
    vector<double> wma12 = decryptVec(wma12_encrypted_sliced, scheme, secretKey);
    vector<double> wma26 = decryptVec(wma26_encrypted, scheme, secretKey);
    vector<double> wma_diff = decryptVec(wma_diff_sliced, scheme, secretKey);
    vector<double> wma9 = decryptVec(wma9_encrypted, scheme, secretKey);
    vector<double> macd = decryptVec(macd_encrypted, scheme, secretKey);
    vector<double> decisions = decryptVec(decisions_encrypted, scheme, secretKey);

    exportData(wma12, "../data/wma12_heaan.csv");
    exportData(wma26, "../data/wma26_heaan.csv");
    exportData(wma_diff, "../data/wma_diff_heaan.csv");
    exportData(wma9, "../data/wma9_heaan.csv");
    exportData(macd, "../data/macd_heaan.csv");
    exportData(decisions, "../data/decisions_heaan.csv");

	cout << "Data Export Finished" << endl;
    cout << "\n";

	return 0;
}
