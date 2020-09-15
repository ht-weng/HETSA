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
// MACD analysis functions
// ***************************************************************************************************

// Return the weighted moving averages of the input Ciphertexts in a moving window
inline Ciphertext ema(vector<Ciphertext>& data, Ciphertext prev_ema, long m, Scheme &scheme, long logq, long logp) {

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
inline Ciphertext getSample(long time, long logq, long logp, Scheme& scheme) {
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
inline void assembleSample(Ciphertext& sample, vector<Ciphertext>& past_prices) {
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

    // Every time we receive the encrypted price, we store it in a vector
    vector<Ciphertext> prices_encrypted;

    // EMA window size
    long win = 12;
    vector<Ciphertext> ema12_encrypted;
    vector<Ciphertext> prev_ema12;
    
    Ciphertext zero_encrypted;
    complex<double> zero;
    zero.real(0.0);
    scheme.encryptSingle(zero_encrypted, zero, logp, logq);

    cout << "EMA Analysis Begins" << endl;
    auto ema_start = high_resolution_clock::now();

    for (int i = 0; i < time_max; i++) {
        // Import encrypted data
        Ciphertext sample = getSample(i, logq, logp, scheme);
        assembleSample(sample, prices_encrypted);

        if (i > win-1) {
            // When the number of data samples is larger than the window size, do EMA analysis recursively
            Ciphertext prev_ema = prev_ema12.back();
            if (prev_ema.logq < (logq_boots+logp)) {
                // When the previous EMA's modulus coefficient logq is too small, do bootstrapping
                scheme.modDownToAndEqual(prev_ema, logq_boots);
                scheme.bootstrapAndEqual(prev_ema, logq_boots, logQ, logT, logI);
                prev_ema.logp = logp;
            }
            Ciphertext new_ema = ema(prices_encrypted, prev_ema, win, scheme, logq, logp);
            prev_ema12.push_back(new_ema);
            // Make sure the encrypted EMAs have the same mudulus coefficient logq
            scheme.modDownToAndEqual(new_ema, logq_boots);
            new_ema.logp = logp;
            ema12_encrypted.push_back(new_ema);
        }
        
        if (i == win-1) {
            // When the number of data samples equals to the window size, calculate the first EMA
            Ciphertext new_ema = ema(prices_encrypted, zero_encrypted, win, scheme, logq, logp);
            prev_ema12.push_back(new_ema);
            // Make sure the encrypted EMAs have the same mudulus coefficient logq
            scheme.modDownToAndEqual(new_ema, logq_boots);
            new_ema.logp = logp;
            ema12_encrypted.push_back(new_ema);
        }
    }

    auto ema_stop = high_resolution_clock::now();
    auto ema_duration = duration_cast<seconds>(ema_stop - ema_start);
	cout << "EMA Analysis Finished" << endl;
    cout << "EMA Analysis Duration:  " << ema_duration.count() << " seconds" << endl;
    cout << "\n";

    cout << "Data Export Begins" << endl;
    auto export_start = high_resolution_clock::now();
 
    vector<double> ema12 = decryptVec(ema12_encrypted, scheme, secretKey);
    exportData(ema12, "../data/ema12_heaan.csv");
    
    auto export_stop = high_resolution_clock::now();
    auto export_duration = duration_cast<seconds>(export_stop - export_start);
	cout << "Data Export Finished" << endl;
    cout << "Data Export Duration:  " << export_duration.count() << " seconds" << endl;
    cout << "\n";

	return 0;
}
