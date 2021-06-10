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
inline Ciphertext getSampleSingle(string file, long time, long logq, long logp, Scheme& scheme) {
    // Read data from csv file
    vector<double> inputs;
	inputs = csv2vec(file);

    complex<double> sample;
    if (time <= inputs.size()) {
        sample.real(inputs[time]);
    } else {
        cout << "Error: time exceeds sample data range! " << endl;
        sample.real(0.0);
    }
    Ciphertext sample_encrypted;
    scheme.encryptSingle(sample_encrypted, sample, logp, logq);
    return sample_encrypted;
}

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

// Calculate the factorial of a number
inline long factorial(int x) {
    long res = 1;
    for (int i = 1; i <= x; ++i) {
        res *= i;
    }
    return res;
}

// inline Ciphertext s0_sin(Ciphertext x, int r, int terms, Scheme &scheme, SchemeAlgo &algo, long logq, long logp, long logn) {
//     long n = pow(2, logn);
//     Ciphertext result;
//     // // Initialise result with 0
//     // complex<double>* zero = new complex<double>[n];
//     // for (int i = 0; i < n; i++) {
//     //     zero[i].real(0.0);
//     // }
//     // scheme.encrypt(result, zero, n, logp, logq);

//     for (int i = 0; i < terms; i++) {
//         Ciphertext new_term;
//         complex<double> coeff;
//         coeff.real((pow(-1, i)/factorial(2*i+1)) * pow((1/pow(2, r)), (2*i+1)));

//         cout << "sin fac: " << factorial(2*i+1) << endl;
//         cout << "sin coeff: " << (pow(-1, i)/factorial(2*i+1)) * pow((1/pow(2, r)), (2*i+1)) << endl;
        
//         algo.power(new_term, x, logp, long(2*i+1));
//         scheme.multByConstAndEqual(new_term, coeff, logp);
//         scheme.reScaleByAndEqual(new_term, logp);

//         // scheme.modDownToAndEqual(result, new_term.logq);
//         if (i == 0) {
//             result = new_term;
//         } else {
//             scheme.modDownToAndEqual(result, new_term.logq);
//             scheme.addAndEqual(result, new_term);
//         }
//         // scheme.add(result, result, new_term);

//         cout << "sin order: " << 2*i+1 << endl;
//     }
//     cout << "return" << endl;
//     return result;
// }


inline vector<Ciphertext> s0_sin(Ciphertext x, int r, int terms, Scheme &scheme, SchemeAlgo &algo, long logq, long logp, long logn) {
    vector<Ciphertext> result_vec;
    Ciphertext result, x_encrypted, x3_encrypted, x5_encrypted;

    // Set up coefficients
    complex<double> coeff0, coeff1, coeff2;
    coeff0.real(0.125);
    coeff1.real(-0.000325520833);
    coeff2.real(0.0000002543131504167);

    // Calculate 0.125*x
    scheme.multByConst(x_encrypted, x, coeff0, logp);
    scheme.reScaleByAndEqual(x_encrypted, logp);

    // Calculate -0.000325520833*x^3
    algo.power(x3_encrypted, x, logp, 3);
    scheme.multByConstAndEqual(x3_encrypted, coeff1, logp);
    scheme.reScaleByAndEqual(x3_encrypted, logp);

    // Calculate 0.0000002543131504167*x^5
    algo.power(x5_encrypted, x, logp, 5);
    scheme.multByConstAndEqual(x5_encrypted, coeff2, logp);
    scheme.reScaleByAndEqual(x5_encrypted, logp);

    // Calculate s0(x)
    scheme.modDownToAndEqual(x3_encrypted, x5_encrypted.logq);
    scheme.add(result, x3_encrypted, x5_encrypted);
    scheme.modDownToAndEqual(x_encrypted, result.logq);
    scheme.addAndEqual(result, x_encrypted);

    result_vec.push_back(result);
    return result_vec;
}


// inline Ciphertext c0(Ciphertext x, int r, int terms, Scheme &scheme, SchemeAlgo &algo, long logq, long logp, long logn) {
//     long n = pow(2, logn);
//     Ciphertext result;
//     // // Initialise result with 0
//     // complex<double>* zero = new complex<double>[n];
//     // for (int i = 0; i < n; i++) {
//     //     zero[i].real(0.0);
//     // }
//     // scheme.encrypt(result, zero, n, logp, logq);
    
//     for (int i = 0; i < terms; i++) {
//         Ciphertext new_term;
//         complex<double> coeff;
//         coeff.real((pow(-1, i)/factorial(2*i)) * pow((1/pow(2, r)), (2*i)));

//         cout << "cos fac: " << factorial(2*i+1) << endl;
//         cout << "cos coeff: " << (pow(-1, i)/factorial(2*i+1)) * pow((1/pow(2, r)), (2*i+1)) << endl;

//         algo.power(new_term, x, logp, long(2*i));
//         scheme.multByConstAndEqual(new_term, coeff, logp);
//         scheme.reScaleByAndEqual(new_term, logp);

//         scheme.modDownToAndEqual(result, new_term.logq);
//         if (i == 0) {
//             result = new_term;
//         } else {
//             scheme.addAndEqual(result, new_term);
//         }
//         // scheme.addAndEqual(result, new_term);

//         cout << "cos order: " << 2*i+1 << endl;
//     }

//     return result;
// }


inline vector<Ciphertext> c0(Ciphertext x, int r, int terms, Scheme &scheme, SchemeAlgo &algo, long logq, long logp, long logn) {
    vector<Ciphertext> result_vec;
    Ciphertext result, x2_encrypted, x4_encrypted;

    // Set up coefficients
    complex<double> coeff0, coeff1, coeff2;
    coeff0.real(0.125);
    coeff1.real(-0.000325520833);
    coeff2.real(0.0000002543131504167);

    // Calculate -0.0078125*x^2
    algo.power(x2_encrypted, x, logp, 2);
    scheme.multByConstAndEqual(x2_encrypted, coeff1, logp);
    scheme.reScaleByAndEqual(x2_encrypted, logp);

    // Calculate 1.017252604167*x^4
    algo.power(x4_encrypted, x, logp, 4);
    scheme.multByConstAndEqual(x4_encrypted, coeff2, logp);
    scheme.reScaleByAndEqual(x4_encrypted, logp);

    // Calculate c0(x)
    scheme.modDownToAndEqual(x2_encrypted, x4_encrypted.logq);
    scheme.add(result, x2_encrypted, x4_encrypted);
    scheme.addConstAndEqual(result, coeff0, logp);

    result_vec.push_back(result);
    return result_vec;
}

inline Ciphertext sin_approx(Ciphertext x, int r, int terms, Scheme &scheme, SchemeAlgo &algo, long logq, 
    long logp, long logn, long logQ, long logT, long logI, long logq_boots) {
    Ciphertext sk, ck;
    vector<Ciphertext> sk_vec, ck_vec;
    complex<double> coeff;

    cout << "Start Calculating s0 and c0" << endl;

    coeff.real(2.0);
    sk_vec = s0_sin(x, r, terms, scheme, algo, logq, logp, logn);
    sk = sk_vec[0];
    ck_vec = c0(x, r, terms, scheme, algo, logq, logp, logn);
    ck = ck_vec[0];

    cout << "Finish Calculating s0 and c0" << endl;

    cout << "Start Initial Bootstrapping" << endl;

    scheme.bootstrapAndEqual(sk, logq_boots, logQ, logT, logI);
    scheme.bootstrapAndEqual(ck, logq_boots, logQ, logT, logI);

    cout << "Finish Initial Bootstrapping" << endl;

    for (int i = 0; i < r; i ++) {
        Ciphertext sk1, ck1, sk_square, ck_square;
        // Calculate ck1
        scheme.square(sk_square, sk);
        scheme.reScaleByAndEqual(sk_square, logp);
        scheme.square(ck_square, ck);
        scheme.reScaleByAndEqual(ck_square, logp);
        scheme.modDownToAndEqual(ck_square, sk_square.logq);
        scheme.modDownToAndEqual(sk_square, ck_square.logq);
        scheme.sub(ck1, ck_square, sk_square);

        // Calculate sk1
        scheme.mult(sk1, sk, ck);
        scheme.reScaleByAndEqual(sk1, logp);
        scheme.multByConstAndEqual(sk1, coeff, logp);
        scheme.reScaleByAndEqual(sk1, logp);

        // Bootstrapping
        cout << "Start " << i << " Bootstrapping" << endl;
        scheme.bootstrapAndEqual(sk1, logq_boots, logQ, logT, logI);
        scheme.bootstrapAndEqual(ck1, logq_boots, logQ, logT, logI);
        cout << "Finish " << i << " Bootstrapping" << endl;

        // Reassign
        sk = sk1;
        ck = ck1;
    }

    return sk;
}


//********************************************************************************
// Main function
//********************************************************************************
int main() {

	// long logq = 600; // Ciphertext polynomial modulus
	// long logp = 30; // Scale
	// long logn = 5; // n = number of slots in a Ciphertext
    long logq = 800; // Ciphertext polynomial modulus
	long logp = 30; // Scale
	long logn = 4; // n = number of slots in a Ciphertext
    // long logq = 300; // Ciphertext polynomial modulus
	// long logp = 25; // Scale
	// long logn = 4; // n = number of slots in a Ciphertext

    long logT = 2;
    long logI = 4;
    long logQ = 1200;
    long logq_boots = 50;

    int r = 3;
    int terms = 3;
    int time_max = 158;

	SetNumThreads(8);
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
    scheme.addBootKey(secretKey, 0, logq_boots + logI);
    SchemeAlgo algo(scheme);

    cout << "Data Import Begins" << endl;
    auto import_start = high_resolution_clock::now();

    vector<Ciphertext> x_samples;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSampleSingle("../../data/pi_samples.csv", i, logq, logp, scheme);
        assembleSample(sample, x_samples);   
    }

    auto import_stop = high_resolution_clock::now();
    auto import_duration = duration_cast<seconds>(import_stop - import_start); 

    cout << "Data Import Finished" << endl;
    cout << "Data Import Duration:  " << import_duration.count() << " seconds" << endl;
	cout << "\n";

    cout << "Sine Interpolation Test Begins" << endl;
    auto sin_start = high_resolution_clock::now();

    vector<Ciphertext> sin_encrypted;
    
    for (int i = 0; i < x_samples.size(); i++) {

        // Ciphertext sin;
        // sin = sin_approx(x_samples[i], r, terms, scheme, algo, logq, logp, logn, logQ, logT, logI, logq_boots);
        // sin = s0_sin(x_samples[i], r, terms, scheme, algo, logq, logp, logn);
        Ciphertext sk, ck, x;
        vector<Ciphertext> sk_vec, ck_vec;
        complex<double> coeff;
        x = x_samples[i];

        cout << "Start Calculating s0 and c0" << endl;

        coeff.real(2.0);
        sk_vec = s0_sin(x, r, terms, scheme, algo, logq, logp, logn);
        sk = sk_vec[0];
        ck_vec = c0(x, r, terms, scheme, algo, logq, logp, logn);
        ck = ck_vec[0];

        // complex<double>* sk_val;
        // sk_val = scheme.decrypt(secretKey, sk);
        // cout << "sk: " << real(sk_val[0]) << endl;

        // complex<double>* ck_val;
        // ck_val = scheme.decrypt(secretKey, ck);
        // cout << "ck: " << real(ck_val[0]) << endl;

        cout << "Finish Calculating s0 and c0" << endl;

        cout << "Start Initial Bootstrapping" << endl;

        scheme.bootstrapAndEqual(sk, logq_boots, logQ, logT, logI);
        cout << "Finish Bootstrapping sk" << endl;
        scheme.bootstrapAndEqual(ck, logq_boots, logQ, logT, logI);
        sk.logp = logp;
        ck.logp = logp;

        cout << "Finish Initial Bootstrapping" << endl;

        // Decrypt
        complex<double>* sk_val;
        sk_val = scheme.decrypt(secretKey, sk);
        cout << "sk: " << real(sk_val[0]) << endl;

        // for (int i = 0; i < r; i ++) {
        //     Ciphertext sk1, ck1, sk_square, ck_square;
        //     // Calculate ck1
        //     scheme.square(sk_square, sk);
        //     scheme.reScaleByAndEqual(sk_square, logp);
        //     scheme.square(ck_square, ck);
        //     scheme.reScaleByAndEqual(ck_square, logp);
        //     scheme.modDownToAndEqual(ck_square, sk_square.logq);
        //     scheme.modDownToAndEqual(sk_square, ck_square.logq);
        //     scheme.sub(ck1, ck_square, sk_square);

        //     // Calculate sk1
        //     scheme.mult(sk1, sk, ck);
        //     scheme.reScaleByAndEqual(sk1, logp);
        //     scheme.multByConstAndEqual(sk1, coeff, logp);
        //     scheme.reScaleByAndEqual(sk1, logp);

        //     // Bootstrapping
        //     cout << "Start " << i << " Bootstrapping" << endl;
        //     scheme.bootstrapAndEqual(sk1, logq_boots, logQ, logT, logI);
        //     scheme.bootstrapAndEqual(ck1, logq_boots, logQ, logT, logI);
        //     sk1.logp = logp;
        //     ck1.logp = logp;
        //     cout << "Finish " << i << " Bootstrapping" << endl;

        //     // Reassign
        //     sk = sk1;
        //     ck = ck1;

        //     cout << "Finish " << i << "Iteration" << endl;
        // }
        // sin_encrypted.push_back(sk);

        cout << "Start 1st Iteration" << endl;
        Ciphertext sk1, ck1, sk_square, ck_square;
        // Calculate ck1
        scheme.square(sk_square, sk);
        scheme.reScaleByAndEqual(sk_square, logp);
        scheme.square(ck_square, ck);
        scheme.reScaleByAndEqual(ck_square, logp);
        scheme.modDownToAndEqual(ck_square, sk_square.logq);
        scheme.modDownToAndEqual(sk_square, ck_square.logq);
        scheme.sub(ck1, ck_square, sk_square);

        // Calculate sk1
        scheme.mult(sk1, sk, ck);
        scheme.reScaleByAndEqual(sk1, logp);
        scheme.multByConstAndEqual(sk1, coeff, logp);
        scheme.reScaleByAndEqual(sk1, logp);

        // Bootstrapping
        cout << "Start 1st Bootstrapping" << endl;
        scheme.bootstrapAndEqual(sk1, logq_boots, logQ, logT, logI);
        scheme.bootstrapAndEqual(ck1, logq_boots, logQ, logT, logI);
        sk1.logp = logp;
        ck1.logp = logp;
        cout << "Finish 1st Bootstrapping" << endl;
        cout << "Finish 1st Iteration" << endl;

        // Decrypt
        complex<double>* sk1_val;
        sk1_val = scheme.decrypt(secretKey, sk1);
        cout << "sk1: " << real(sk1_val[0]) << endl;


        cout << "Start 2nd Iteration" << endl;
        Ciphertext sk2, ck2, sk_square2, ck_square2;
        // Calculate ck2
        scheme.square(sk_square2, sk1);
        scheme.reScaleByAndEqual(sk_square2, logp);
        scheme.square(ck_square2, ck1);
        scheme.reScaleByAndEqual(ck_square2, logp);
        scheme.modDownToAndEqual(ck_square2, sk_square2.logq);
        scheme.modDownToAndEqual(sk_square2, ck_square2.logq);
        scheme.sub(ck2, ck_square2, sk_square2);

        // Calculate sk2
        scheme.mult(sk2, sk1, ck1);
        scheme.reScaleByAndEqual(sk2, logp);
        scheme.multByConstAndEqual(sk2, coeff, logp);
        scheme.reScaleByAndEqual(sk2, logp);

        // Bootstrapping
        cout << "Start 2nd Bootstrapping" << endl;
        scheme.bootstrapAndEqual(sk2, logq_boots, logQ, logT, logI);
        scheme.bootstrapAndEqual(ck2, logq_boots, logQ, logT, logI);
        sk2.logp = logp;
        ck2.logp = logp;
        cout << "Finish 2nd Bootstrapping" << endl;
        cout << "Finish 2nd Iteration" << endl;

        // Decrypt
        complex<double>* sk2_val;
        sk2_val = scheme.decrypt(secretKey, sk2);
        cout << "sk2: " << real(sk2_val[0]) << endl;


        cout << "Start 3rd Iteration" << endl;
        Ciphertext sk3, ck3, sk_square3, ck_square3;
        // Calculate ck3
        scheme.square(sk_square3, sk2);
        scheme.reScaleByAndEqual(sk_square3, logp);
        scheme.square(ck_square3, ck2);
        scheme.reScaleByAndEqual(ck_square3, logp);
        scheme.modDownToAndEqual(ck_square3, sk_square3.logq);
        scheme.modDownToAndEqual(sk_square3, ck_square3.logq);
        scheme.sub(ck3, ck_square3, sk_square3);

        // Calculate sk3
        scheme.mult(sk3, sk2, ck2);
        scheme.reScaleByAndEqual(sk3, logp);
        scheme.multByConstAndEqual(sk3, coeff, logp);
        scheme.reScaleByAndEqual(sk3, logp);

        // Bootstrapping
        cout << "Start 3rd Bootstrapping" << endl;
        scheme.bootstrapAndEqual(sk3, logq_boots, logQ, logT, logI);
        scheme.bootstrapAndEqual(ck3, logq_boots, logQ, logT, logI);
        sk3.logp = logp;
        ck3.logp = logp;
        cout << "Finish 3rd Bootstrapping" << endl;
        cout << "Finish 3rd Iteration" << endl;

        // Decrypt
        complex<double>* sk3_val;
        sk3_val = scheme.decrypt(secretKey, sk3);
        cout << "sk3: " << real(sk3_val[0]) << endl;

        // sin_encrypted.push_back(sk3);
    }

    auto sin_stop = high_resolution_clock::now();
    auto sin_duration = duration_cast<seconds>(sin_stop - sin_start);
    cout << "Sine Interpolation Test Finished" << endl;
    cout << "Sine Interpolation Test Duration:  " << sin_duration.count() << " seconds" << endl;
    cout << "\n";

    cout << "Data Export Begins" << endl;
    auto export_start = high_resolution_clock::now();

    vector<double> sin_plain = decryptVec(sin_encrypted, scheme, secretKey);
    exportData(sin_plain, "../../data/sine_interpolated_heaan.csv");

    auto export_stop = high_resolution_clock::now();
    auto export_duration = duration_cast<seconds>(export_stop - export_start);
	cout << "Data Export Finished" << endl;
    cout << "Data Export Duration:  " << export_duration.count() << " seconds" << endl;
    cout << "\n";

	return 0;
}
