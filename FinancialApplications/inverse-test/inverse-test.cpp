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
// Data input/output functions
// ***************************************************************************************************

// Simulate the data import process. Return corresponding sample data according to the time
inline Ciphertext getSample(string file, long time, long logq, long logp, long logn, Scheme& scheme) {
    // Read data from csv file
    vector<double> inputs;
	inputs = csv2vec(file);

    long n = pow(2, logn);
    complex<double>* sample = new complex<double>[n];
    if (time <= inputs.size()) { 
        sample[0].real(inputs[time]);
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
    delete[] sample;
    return sample_encrypted;
}

// Append new sample data to the history data
inline void assembleSample(Ciphertext& sample, vector<Ciphertext>& past_inputs) {
    past_inputs.push_back(sample);
}

// Decrytp a vector of Ciphertexts into a vector of doubles
inline vector<double> decryptVec(vector<Ciphertext>& ct, Scheme& scheme, SecretKey& secretKey) {
    vector<double> res;
    for (int i = 0; i < ct.size(); i++) {
        complex<double>* val;
        val = scheme.decrypt(secretKey, ct[i]);
        res.push_back(real(val[0]));
        delete val;
    }
    return res;
}

// Write data to csv files for plotting
inline void exportData(vector<double> vec, string file) {
    ofstream output_file(file);
    ostream_iterator<double> output_iterator(output_file, "\n");
    copy(vec.begin(), vec.end(), output_iterator);
}

inline vector<Ciphertext> inverse_series(vector<Ciphertext>& data, Scheme &scheme, SchemeAlgo &algo, long logq, long logp, long steps) {
    vector<Ciphertext> result;

    for (int i = 0; i < data.size(); i++) {
        Ciphertext inv;
        algo.inverse(inv, data[i], logp, steps);
        result.push_back(inv);
    }

    return result;
}

//********************************************************************************
// Main function
//********************************************************************************
int main() {

	long logq = 300; // Ciphertext polynomial modulus
	long logp = 30; // Scale
	long logn = 5; // n = number of slots in a Ciphertext
    // long logq = 300; // Ciphertext polynomial modulus
	// long logp = 25; // Scale
	// long logn = 4; // n = number of slots in a Ciphertext
	long time_max = 200; // Sample data size
    long steps = 6;

	SetNumThreads(8);
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
    SchemeAlgo algo(scheme);

    cout << "Data Import Begins" << endl;
    auto import_start = high_resolution_clock::now();

    // vector<Ciphertext> prices_encrypted;
    // for (int i = 0; i < time_max; i++) {
    //     Ciphertext sample = getSample("../data/apple_prices.csv", i, logq, logp, logn, scheme);
    //     assembleSample(sample, prices_encrypted);   
    // }

    // vector<Ciphertext> noises_encrypted;
    // for (int i = 0; i < time_max; i++) {
    //     Ciphertext sample = getSample("../data/white_noise.csv", i, logq, logp, logn, scheme);
    //     assembleSample(sample, noises_encrypted);   
    // }

    vector<Ciphertext> small_noises_encrypted;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample("../data/white_noise_small.csv", i, logq, logp, logn, scheme);
        assembleSample(sample, small_noises_encrypted);   
    }


    auto import_stop = high_resolution_clock::now();
    auto import_duration = duration_cast<seconds>(import_stop - import_start); 

    cout << "Data Import Finished" << endl;
    cout << "Data Import Duration:  " << import_duration.count() << " seconds" << endl;
	cout << "\n";


    // cout << "Inverse Test on Stock Prices Begins" << endl;
    // auto inv_start_1 = high_resolution_clock::now();

    // vector<Ciphertext> prices_inv = inverse_series(prices_encrypted, scheme, algo, logq, logp, steps);

    // auto inv_stop_1 = high_resolution_clock::now();
    // auto inv_duration_1 = duration_cast<seconds>(inv_stop_1 - inv_start_1);
    // cout << "Inverse Test on Stock Prices Finished" << endl;
    // cout << "Inverse Test on Stock Prices Duration:  " << inv_duration_1.count() << " seconds" << endl;
    // cout << "\n";

    // cout << "Inverse Test on White Noises Begins" << endl;
    // auto inv_start_2 = high_resolution_clock::now();

    // vector<Ciphertext> noises_inv = inverse_series(noises_encrypted, scheme, algo, logq, logp, steps);

    // auto inv_stop_2 = high_resolution_clock::now();
    // auto inv_duration_2 = duration_cast<seconds>(inv_stop_2 - inv_start_2);
    // cout << "Inverse Test on White Noises Finished" << endl;
    // cout << "Inverse Test on White Noises Duration:  " << inv_duration_2.count() << " seconds" << endl;
    // cout << "\n";


    cout << "Inverse Test on Small White Noises Begins" << endl;
    auto inv_start_3 = high_resolution_clock::now();

    vector<Ciphertext> small_noises_inv = inverse_series(small_noises_encrypted, scheme, algo, logq, logp, steps);

    auto inv_stop_3 = high_resolution_clock::now();
    auto inv_duration_3 = duration_cast<seconds>(inv_stop_3 - inv_start_3);
    cout << "Inverse Test on Small White Noises Finished" << endl;
    cout << "Inverse Test on Small White Noises Duration:  " << inv_duration_3.count() << " seconds" << endl;
    cout << "\n";

    cout << "Data Export Begins" << endl;
    auto export_start = high_resolution_clock::now();
    
    // vector<double> prices_inv_plain = decryptVec(prices_inv, scheme, secretKey);
    // exportData(prices_inv_plain, "../data/prices_inv_plain.csv");

    // vector<double> noises_inv_plain = decryptVec(noises_inv, scheme, secretKey);
    // exportData(noises_inv_plain, "../data/noises_inv_plain.csv");

    vector<double> noises_inv_plain = decryptVec(small_noises_inv, scheme, secretKey);
    exportData(noises_inv_plain, "../data/small_noises_inv_plain.csv");

    auto export_stop = high_resolution_clock::now();
    auto export_duration = duration_cast<seconds>(export_stop - export_start);
	cout << "Data Export Finished" << endl;
    cout << "Data Export Duration:  " << export_duration.count() << " seconds" << endl;
    cout << "\n";

	return 0;
}
