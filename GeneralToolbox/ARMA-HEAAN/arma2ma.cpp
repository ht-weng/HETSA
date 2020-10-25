/*
* Autoregressive model with HEAAN Library
* The results are output into the data folder
* Run ema-heaan.ipynb to visualise the results
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
#include <chrono>
#include <iterator>
#include "../../HEAAN/src/HEAAN.h"

using namespace std;
using namespace std::chrono;

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
// ARMA to MA
// ***************************************************************************************************

// Convert ARMA model coefficients into MA model coefficients
// thetas: MA coefficients
// phis: AR coefficients
inline vector<Ciphertext> arma2ma(vector<Ciphertext>& thetas, vector<Ciphertext>& phis, int order, Scheme &scheme, long logq, long logp, 
    long logq_boots, long logQ, long logT, long logI) {

    int t = thetas.size();
    int p = phis.size();

    Ciphertext zero_encrypted;
    complex<double> zero;
    zero.real(0.0);
    scheme.encryptSingle(zero_encrypted, zero, logp, logq);
    Ciphertext one_encrypted;
    complex<double> one;
    one.real(1.0);
    scheme.encryptSingle(one_encrypted, one, logp, logq);

    // initialise result with zeros
    vector<Ciphertext> result;
    result.push_back(one_encrypted);
    for (int i = 1; i < order; i++) {
        result.push_back(zero_encrypted);
    }

    // negate elements in phis
    vector<Ciphertext> n_phis;
    for (int i = 0; i < p; i++) {
        Ciphertext n_phis_coeff;
        scheme.sub(n_phis_coeff, zero_encrypted, phis[i]);
        n_phis.push_back(n_phis_coeff);
    }
    
    // padding
    if (t < order) {
        for (int i = t; i < order; i++) {
            thetas.push_back(zero_encrypted);
        }
    }
    if (p < order) {
        for (int i = p; i < order; i++) {
            n_phis.push_back(zero_encrypted);
        }
    }
    
    for (int i = 0; i < order-1; i++) {
        result[i+1] = thetas[i+1];

        for (int j = 0; j < (min(i, p)+1); j++) {
            // issue: multiplication levels of result[i+1] and result[i-j]*n_phis[j+1] must be the same
            // to allow addition
            // calculating result is a recursive process

            // result[i+1] = result[i+1] + result[i-j]*n_phis[j+1];

            Ciphertext tmp;
            scheme.mult(tmp, result[i+1], n_phis[j+1]);
            scheme.addAndEqual(result[i+1], tmp);
        }
    }
    return result;
    
}

// Return the weighted moving averages of the input Ciphertexts in a moving window
Ciphertext ema(vector<Ciphertext>& data, Ciphertext prev_ema, long m, Scheme &scheme, long logq, long logp) {

    complex<double> multiplier1;
    multiplier1.real(2.0/float(m+1));
    complex<double> multiplier2;
    multiplier2.real(1.0-multiplier1.real());
    complex<double> denom;
    denom.real(1.0/float(m));

    if (data.size() < m) {
        Ciphertext zero_encrypted;
        complex<double> zero;
        zero.real(0.0);
        scheme.encryptSingle(zero_encrypted, zero, logp, logq);
        cout << "Error: the length of stock prices is smaller than the window size!" << endl;
        return zero_encrypted;
    } else if (data.size() == m) {
        // The first EMA is the average of the first m stock prices
        Ciphertext new_ema = sum(data, scheme);
        scheme.multByConstAndEqual(new_ema, denom, logp);
        scheme.reScaleByAndEqual(new_ema, logp);
        return new_ema;
    } else {
        // EMA(current) = ((Price(current) - EMA(prev)) * Multiplier) + EMA(prev)
        // = (Price(current) * Multiplier + EMA(prev) * (1 - Multiplier)
        Ciphertext curr_price = data.back();
        scheme.multByConstAndEqual(curr_price, multiplier1, logp);
        scheme.reScaleByAndEqual(curr_price, logp);
        scheme.multByConstAndEqual(prev_ema, multiplier2, logp);
        scheme.reScaleByAndEqual(prev_ema, logp);
        if (curr_price.logq > prev_ema.logq) {
            scheme.modDownToAndEqual(curr_price, prev_ema.logq);
        } else if (curr_price.logq < prev_ema.logq) {
            scheme.modDownToAndEqual(prev_ema, curr_price.logq);
        }
        Ciphertext new_ema;
        scheme.add(new_ema, curr_price, prev_ema);
        return new_ema;
    }
    
}

// ***************************************************************************************************
// Data input/output functions
// ***************************************************************************************************

// Simulate the data import process. Return corresponding sample data according to the time
Ciphertext getSample(long time, long logq, long logp, Scheme& scheme) {
    // Read data from csv file
    vector<double> prices;
	prices = csv2vec("../data/apple_prices.csv");

    complex<double> sample;
    if (time <= prices.size()) {
        sample.real(prices[time]);
    } else {
        cout << "Error: time exceeds sample data range! " << endl;
        sample.real(0.0);
    }
    Ciphertext sample_encrypted;
    scheme.encryptSingle(sample_encrypted, sample, logp, logq);
    return sample_encrypted;
}

// Append new sample data to the history data
void assembleSample(Ciphertext& sample, vector<Ciphertext>& past_prices) {
    past_prices.push_back(sample);
}

// Decrytp a vector of Ciphertexts into a vector of doubles
inline vector<double> decryptVec(vector<Ciphertext>& ct, Scheme& scheme, SecretKey& secretKey) {
    vector<double> res;
    for (int i = 0; i < ct.size(); i++) {
        complex<double> val = scheme.decryptSingle(secretKey, ct[i]);
        res.push_back(real(val));
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
	// long logn = 5; // n = number of slots in a Ciphertext
    long logT = 2;
    long logI = 4;
    long logQ = 1200;
    long logq_boots = 50;

    long time_max = 100; // Sample data size

	SetNumThreads(8);
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
    scheme.addBootKey(secretKey, 0, logq_boots + logI);

    // set ARMA parameters and order for ARMA2MA
    vector<double> maparams{0.35};
    vector<double> arparams{1.0, -0.65};
    int order = 10;

    // thetas: MA coefficients
    // phis: AR coefficients
    vector<Ciphertext> thetas;
    vector<Ciphertext> phis;

    // Encrypt ARMA coefficients
    for (int i = 0; i < maparams.size(); i++) {
        Ciphertext coeff_encrypted;
        scheme.encryptSingle(coeff_encrypted, maparams[i], logp, logq);
        thetas.push_back(coeff_encrypted);
    }

    for (int i = 0; i < arparams.size(); i++) {
        Ciphertext coeff_encrypted;
        scheme.encryptSingle(coeff_encrypted, arparams[i], logp, logq);
        phis.push_back(coeff_encrypted);
    }

    vector<Ciphertext> result_encrypted;
    auto arma2ma_start = high_resolution_clock::now();

    // ARMA2MA
    result_encrypted = arma2ma(thetas, phis, order, scheme, logq, logp, logq_boots, logQ, logT, logI);

    auto arma2ma_stop = high_resolution_clock::now();
    auto arma2ma_duration = duration_cast<seconds>(arma2ma_stop - arma2ma_start);
    cout << "ARMA2MA Duration:  " << arma2ma_duration.count() << " seconds" << endl;
    cout << endl;

    // Decrypt and print resulting coefficients
    vector<double> result = decryptVec(result_encrypted, scheme, secretKey);
    printVector(result);

	return 0;
}
