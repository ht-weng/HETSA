/*
* Moving Average Convergence Divergence with SEAL Library
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
#include <chrono>
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

//********************************************************************************
// SEAL helper functions
//********************************************************************************
/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(std::shared_ptr<seal::SEALContext> context)
{
    // Verify parameters
    if (!context)
    {
        throw std::invalid_argument("context is not set");
    }
    auto &context_data = *context->key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::BFV:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::CKKS:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " <<
        context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_mod_count = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_mod_count - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::BFV)
    {
        std::cout << "|   plain_modulus: " << context_data.
            parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
inline std::ostream &operator <<(std::ostream &stream, seal::parms_id_type parms_id)
{
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    stream << std::hex << std::setfill('0')
        << std::setw(16) << parms_id[0] << " "
        << std::setw(16) << parms_id[1] << " "
        << std::setw(16) << parms_id[2] << " "
        << std::setw(16) << parms_id[3] << " ";

    std::cout.copyfmt(old_fmt);

    return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template<typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if(slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if(vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    std::cout.copyfmt(old_fmt);
}

//********************************************************************************
// MACD analysis functions
//********************************************************************************

// Return the sum of all the Ciphertexts in a vector
inline Ciphertext sum(vector<Ciphertext>& vec, CKKSEncoder &encoder, Encryptor &encryptor, 
    Evaluator &evaluator, double scale) {
    Ciphertext result_encrypted = vec[0];

    for(int i = 1; i < vec.size(); i++) {
        Ciphertext val = vec[i];
        evaluator.add_inplace(result_encrypted, val);
    }
    return result_encrypted;
}

// Return a subset of the vector of Ciphertexts
inline vector<Ciphertext> slice(vector<Ciphertext>& vec, int start=0, int end=-1) {
    int oldlen = vec.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }
    vector<Ciphertext> res;

    for (int i = 0; i < newlen; i++) {
        res.push_back(vec[start+i]);
    }
    return res;
}

// Import data in csv file into a vector of doubles
inline vector<double> csv2vec(string file) {
    vector<double> data;
    ifstream inputFile(file);
    int l = 0;

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

// Return the weighted moving averages of the input Ciphertexts in a moving window
inline vector<Ciphertext> wma(vector<Ciphertext>& data, int n, int level, CKKSEncoder &encoder, 
    Encryptor &encryptor, Evaluator &evaluator, double scale, parms_id_type *parms_ids) {
    vector<Ciphertext> wma;
    vector<Plaintext> weights;

    // Generate a list of weights of the window size
    for (int i = 0; i < n; i++) {
        vector<double> w;
        Plaintext w_plain;
        w.push_back(2.0*(i+1.0)/(n*(n+1.0)));
        encoder.encode(w, scale, w_plain);
        weights.push_back(w_plain);
    }

    // Multiply the corresponding Ciphertext and weight in a window and sum the results up 
    for (int i = 0; i < data.size()-n; i++) {
        vector<Ciphertext> data_sliced;
        vector<Ciphertext> window;

        // Get the data in the moving window
        data_sliced = slice(data, i, i+n);

        for (int j = 0; j < n; j++) {
            Ciphertext tmp = data_sliced[j];
            Plaintext tmp_weight = weights[j];

            // Multiply the Ciphertext and weight and rescale the result
            evaluator.mod_switch_to_inplace(tmp_weight, parms_ids[level]);
            evaluator.multiply_plain_inplace(tmp, tmp_weight);
            evaluator.rescale_to_next_inplace(tmp);
            tmp.scale() = scale;
            window.push_back(tmp);
        }

        // Sum the multiplication results up to get the weighted moving average
        Ciphertext res = sum(window, encoder, encryptor, evaluator, scale);
        wma.push_back(res);
    }
    return wma;
}

// ***************************************************************************************************
// Old approximation approach
// ***************************************************************************************************
// // Return the trading decisions based on the MACD signals
// // Denote MACD signals as m(t) where t is the time and sign function as sign()
// // Decision(t) = (m(t-1)-m(t))*(sign(m(t-1)*m(t))-1)
// // sign(m(t-1)*m(t))-1 indicates the turning point, where MACD signals cross the X axis
// // m(t-1)-m(t) indicates the trend
// // Decision(t) = 0 means "Do nothing" or "Hold"
// // Decision(t) > 0 means "Buy"
// // Decision(t) < 0 means "Sell"
// // We use a polynomial approximation of tanh to approximate sign function
// // The approximation of sign function: f(x) = 0.375*x^5 - 1.25*x^3 + 1.875*x where x = m(t)*m(t-1)
// inline vector<Ciphertext> decision(vector<Ciphertext> macd_encrypted, int level, CKKSEncoder &encoder, 
//         Encryptor &encryptor, Evaluator &evaluator, RelinKeys &relin_keys, double scale, parms_id_type *parms_ids) {
    
//     vector<Ciphertext> decisions_encrypted;
//     for (int i = 1; i < macd_encrypted.size(); i++) {

//         // Store the latest two MACD signals m(t) and m(t-1)
//         Ciphertext mt = macd_encrypted[i];
//         Ciphertext mt_1 = macd_encrypted[i-1];

//         // Set up coefficients
//         Plaintext coeff1_plain, coeff2_plain, coeff3_plain, offset_plain, coeff_norm_plain;
//         Ciphertext offset_encrypted;
//         encoder.encode(0.375, scale, coeff1_plain);
//         encoder.encode((-1.25), scale, coeff2_plain);
//         encoder.encode(1.875, scale, coeff3_plain);
//         encoder.encode(-1.0, scale, offset_plain);
//         encryptor.encrypt(offset_plain, offset_encrypted);
//         encoder.encode(0.25, scale, coeff_norm_plain);
//         evaluator.mod_switch_to_inplace(coeff_norm_plain, parms_ids[level]);

//         // Normalise the MACD signals so that most of the data is in range [-1, 1]
//         evaluator.multiply_plain_inplace(mt, coeff_norm_plain);
//         evaluator.rescale_to_next_inplace(mt);
//         mt.scale() = scale;
//         evaluator.multiply_plain_inplace(mt_1, coeff_norm_plain);
//         evaluator.rescale_to_next_inplace(mt_1);
//         mt_1.scale() = scale;

//         // Calculate m(t-1)-m(t) and m(t-1)*m(t)
//         // m(t-1)*m(t) is denoted as x in the following steps
//         Ciphertext recent_diff;
//         evaluator.sub(mt_1, mt, recent_diff);
//         Ciphertext x;
//         evaluator.multiply(mt_1, mt, x);
//         evaluator.relinearize_inplace(x, relin_keys);
//         evaluator.rescale_to_next_inplace(x);
//         x.scale() = scale;

//         // Calculate x^2, which is at level 5
//         Ciphertext x2_encrypted;
//         evaluator.square(x, x2_encrypted);
//         evaluator.relinearize_inplace(x2_encrypted, relin_keys);
//         evaluator.rescale_to_next_inplace(x2_encrypted);
//         x2_encrypted.scale() = scale;

//         // Calculate x^3, which is at level 6
//         // x^2 is at level 5, switch mod to ensure parameters id match
//         Ciphertext x3_encrypted;
//         evaluator.mod_switch_to_inplace(x, parms_ids[level+3]);
//         evaluator.multiply(x, x2_encrypted, x3_encrypted);
//         evaluator.relinearize_inplace(x3_encrypted, relin_keys);
//         evaluator.rescale_to_next_inplace(x3_encrypted);
//         x3_encrypted.scale() = scale;

//         // Calculate x^5, which is at level 7
//         // x^3 is at level 6, switch mod to ensure parameters id match
//         Ciphertext x5_encrypted;        
//         evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[level+4]);
//         evaluator.multiply(x2_encrypted, x3_encrypted, x5_encrypted);
//         evaluator.relinearize_inplace(x5_encrypted, relin_keys);
//         evaluator.rescale_to_next_inplace(x5_encrypted);
//         x5_encrypted.scale() = scale;

//         // Calculate 0.375*x^5, which is at level 8
//         // x^5 is at level 7, switch mod to ensure parameters id match
//         evaluator.mod_switch_to_inplace(coeff1_plain, parms_ids[level+5]);
//         evaluator.multiply_plain_inplace(x5_encrypted, coeff1_plain);
//         evaluator.rescale_to_next_inplace(x5_encrypted);
//         x5_encrypted.scale() = scale;

//         // Calculate -1.25*x^3, which is at level 7
//         // x^3 is at level 6, switch mod to ensure parameters id match
//         evaluator.mod_switch_to_inplace(coeff2_plain, parms_ids[level+4]);
//         evaluator.multiply_plain_inplace(x3_encrypted, coeff2_plain);
//         evaluator.rescale_to_next_inplace(x3_encrypted);
//         x3_encrypted.scale() = scale;

//         // Calculate 1.875*x
//         // Switch mod to ensure parameters id match
//         // Note that previously we haved switched x's level to 5
//         evaluator.mod_switch_to_inplace(coeff3_plain, parms_ids[level+3]);
//         evaluator.multiply_plain_inplace(x, coeff3_plain);
//         evaluator.rescale_to_next_inplace(x);
//         x.scale() = scale;

//         // To do the addition, we have to ensure all the terms have same parms_id and scale
//         evaluator.mod_switch_to_inplace(x3_encrypted, parms_ids[level+6]);
//         evaluator.mod_switch_to_inplace(x, parms_ids[level+6]);
//         evaluator.mod_switch_to_inplace(offset_encrypted, parms_ids[level+6]);

//         // Calculate 0.375*x^5-1.25*x^3+1.875*x-1
//         Ciphertext sign_encrypted;
//         evaluator.add(x5_encrypted, x3_encrypted, sign_encrypted);
//         evaluator.add_inplace(sign_encrypted, x);
//         evaluator.add_inplace(sign_encrypted, offset_encrypted);

//         // Calculate (m(t-1)-m(t))*(f(x)-1)
//         Ciphertext result_encrypted;
//         evaluator.mod_switch_to_inplace(recent_diff, parms_ids[level+6]);
//         sign_encrypted.scale() = scale;
//         recent_diff.scale() = scale;
//         result_encrypted.scale() = scale;
//         evaluator.multiply(sign_encrypted, recent_diff, result_encrypted);
//         evaluator.relinearize_inplace(result_encrypted, relin_keys);
//         evaluator.rescale_to_next_inplace(result_encrypted);
//         result_encrypted.scale() = scale;

//         decisions_encrypted.push_back(result_encrypted);
//     }
//     return decisions_encrypted;
// }

// ***************************************************************************************************
// New approximation approach
// ***************************************************************************************************
// Return the trading decisions based on the MACD signals
// Denote MACD signals as m(t) where t is the time and sign function as sign()
// Decision(t) = (m(t-1)-m(t))*(sign(m(t-1)*m(t))-1)
// sign(m(t-1)*m(t))-1 indicates the turning point, where MACD signals cross the X axis
// m(t-1)-m(t) indicates the trend
// Decision(t) = 0 means "Do nothing" or "Hold"
// Decision(t) > 0 means "Buy"
// Decision(t) < 0 means "Sell"
// We use a polynomial approximation of tanh to approximate sign function
// The approximation of sign function: f(x) = -0.0001*x^9+0.0003*x^8+0.0025*x^7-0.009*x^6 \
    -0.0253*x^5+0.0984*x^4+0.0882*x^3-0.5173*x^2+0.4475*x-0.0753
// where x = m(t)*m(t-1)
inline vector<Ciphertext> decision(vector<Ciphertext> macd_encrypted, int level, double norm, CKKSEncoder &encoder, 
        Encryptor &encryptor, Evaluator &evaluator, RelinKeys &relin_keys, double scale, parms_id_type *parms_ids) {
    
     // Set up coefficients
    Plaintext coeff0_plain, coeff1_plain, coeff2_plain, coeff3_plain, coeff4_plain, coeff5_plain, coeff6_plain,
        coeff7_plain, coeff8_plain, coeff9_plain, coeff_norm_plain;
    encoder.encode((-0.0753), scale, coeff0_plain);
    encoder.encode(0.4475, scale, coeff1_plain);
    encoder.encode((-0.5173), scale, coeff2_plain);
    encoder.encode(0.0882, scale, coeff3_plain);
    encoder.encode(0.0984, scale, coeff4_plain);
    encoder.encode((-0.0253), scale, coeff5_plain);
    encoder.encode((-0.009), scale, coeff6_plain);
    encoder.encode(0.0025, scale, coeff7_plain);
    encoder.encode(0.0003, scale, coeff8_plain);
    encoder.encode((-0.0001), scale, coeff9_plain);
    encoder.encode(norm, scale, coeff_norm_plain);

    vector<Ciphertext> decisions_encrypted;
    
    for (int i = 1; i < macd_encrypted.size(); i++) {

        // Store the latest two MACD signals m(t) and m(t-1)
        Ciphertext mt = macd_encrypted[i];
        Ciphertext mt_1 = macd_encrypted[i-1];

        // Normalise the MACD signals so that most of the data is in range [-1, 1]
        evaluator.mod_switch_to_inplace(coeff_norm_plain, parms_ids[level]);
        evaluator.multiply_plain_inplace(mt, coeff_norm_plain);
        evaluator.rescale_to_next_inplace(mt);
        mt.scale() = scale;
        evaluator.multiply_plain_inplace(mt_1, coeff_norm_plain);
        evaluator.rescale_to_next_inplace(mt_1);
        mt_1.scale() = scale;

        // Calculate m(t-1)-m(t) and m(t-1)*m(t)
        // m(t-1)*m(t) is denoted as x in the following steps, which is at level 4
        Ciphertext recent_diff;
        evaluator.sub(mt_1, mt, recent_diff);
        Ciphertext x_encrypted;
        evaluator.multiply(mt_1, mt, x_encrypted);
        evaluator.relinearize_inplace(x_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x_encrypted);
        x_encrypted.scale() = scale;

        // Calculate x^2, which is at level 5
        Ciphertext x2_encrypted;
        evaluator.square(x_encrypted, x2_encrypted);
        evaluator.relinearize_inplace(x2_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x2_encrypted);
        x2_encrypted.scale() = scale;

        // Calculate x^3, which is at level 6
        // x is at level 4, x^2 is at level 5, switch mod of x to level 5
        Ciphertext x3_encrypted;
        evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[level+3]);
        evaluator.multiply(x_encrypted, x2_encrypted, x3_encrypted);
        evaluator.relinearize_inplace(x3_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x3_encrypted);
        x3_encrypted.scale() = scale;

        // Calculate x^4, which is at level 6
        Ciphertext x4_encrypted;
        evaluator.square(x2_encrypted, x4_encrypted);
        evaluator.relinearize_inplace(x4_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x4_encrypted);
        x4_encrypted.scale() = scale;

        // Calculate x^5, which is at level 7
        // x^2 is at level 5, x^3 is at level 6, swith mod of x^2 to level 6
        Ciphertext x5_encrypted;
        evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[level+4]);
        evaluator.multiply(x3_encrypted, x2_encrypted, x5_encrypted);
        evaluator.relinearize_inplace(x5_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x5_encrypted);
        x5_encrypted.scale() = scale;

        // Calculate x^6, which is at level 7
        // x^4 is at level 6, x^2 is at level 6
        Ciphertext x6_encrypted;
        evaluator.multiply(x2_encrypted, x4_encrypted, x6_encrypted);
        evaluator.relinearize_inplace(x6_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x6_encrypted);
        x6_encrypted.scale() = scale;

        // Calculate x^7, which is at level 7
        // x^4 is at level 6, x^3 is at level 6
        Ciphertext x7_encrypted;
        evaluator.multiply(x3_encrypted, x4_encrypted, x7_encrypted);
        evaluator.relinearize_inplace(x7_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x7_encrypted);
        x7_encrypted.scale() = scale;

        // Calculate x^8, which is at level 7
        Ciphertext x8_encrypted;        
        evaluator.multiply(x4_encrypted, x4_encrypted, x8_encrypted);
        evaluator.relinearize_inplace(x8_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x8_encrypted);
        x8_encrypted.scale() = scale;

        // Calculate x^9, which is at level 8
        // x^5 is at level 7, x^4 is at level 6, switch mod of x^4 to level 7
        Ciphertext x9_encrypted;
        evaluator.mod_switch_to_inplace(x4_encrypted, parms_ids[level+5]);
        evaluator.multiply(x4_encrypted, x5_encrypted, x9_encrypted);
        evaluator.relinearize_inplace(x9_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x9_encrypted);
        x9_encrypted.scale() = scale;

        // Calculate -0.0001*x^9, which is at level 9
        // x^9 is at level 8
        evaluator.mod_switch_to_inplace(coeff9_plain, parms_ids[level+6]);
        evaluator.multiply_plain_inplace(x9_encrypted, coeff9_plain);
        evaluator.rescale_to_next_inplace(x9_encrypted);
        x9_encrypted.scale() = scale;

        // Calculate 0.0003*x^8, which is at level 8
        // x^8 is at level 7
        evaluator.mod_switch_to_inplace(coeff8_plain, parms_ids[level+5]);
        evaluator.multiply_plain_inplace(x8_encrypted, coeff8_plain);
        evaluator.rescale_to_next_inplace(x8_encrypted);
        x8_encrypted.scale() = scale;

        // Calculate 0.0025*x^7, which is at level 8
        // x^7 is at level 7
        evaluator.mod_switch_to_inplace(coeff7_plain, parms_ids[level+5]);
        evaluator.multiply_plain_inplace(x7_encrypted, coeff7_plain);
        evaluator.rescale_to_next_inplace(x7_encrypted);
        x7_encrypted.scale() = scale;

        // Calculate -0.009*x^6, which is at level 8
        // x^6 is at level 7
        evaluator.mod_switch_to_inplace(coeff6_plain, parms_ids[level+5]);
        evaluator.multiply_plain_inplace(x6_encrypted, coeff6_plain);
        evaluator.rescale_to_next_inplace(x6_encrypted);
        x6_encrypted.scale() = scale;

        // Calculate -0.0253*x^5, which is at level 8
        // x^5 is at level 7
        evaluator.mod_switch_to_inplace(coeff5_plain, parms_ids[level+5]);
        evaluator.multiply_plain_inplace(x5_encrypted, coeff5_plain);
        evaluator.rescale_to_next_inplace(x5_encrypted);
        x5_encrypted.scale() = scale;

        // Calculate 0.0984*x^4, which is at level 8
        // x^4 is at level 7 (we switched the mod of x^4 to level 7)
        evaluator.mod_switch_to_inplace(coeff4_plain, parms_ids[level+5]);
        evaluator.multiply_plain_inplace(x4_encrypted, coeff4_plain);
        evaluator.rescale_to_next_inplace(x4_encrypted);
        x4_encrypted.scale() = scale;

        // Calculate 0.0882*x^3, which is at level 7
        // x^3 is at level 6
        evaluator.mod_switch_to_inplace(coeff3_plain, parms_ids[level+4]);
        evaluator.multiply_plain_inplace(x3_encrypted, coeff3_plain);
        evaluator.rescale_to_next_inplace(x3_encrypted);
        x3_encrypted.scale() = scale;

        // Calculate -0.5173*x^2, which is at level 7
        // x^2 is at level 6 (we switched the mod of x^2 to level 6)
        evaluator.mod_switch_to_inplace(coeff2_plain, parms_ids[level+4]);
        evaluator.multiply_plain_inplace(x2_encrypted, coeff2_plain);
        evaluator.rescale_to_next_inplace(x2_encrypted);
        x2_encrypted.scale() = scale;

        // Calculate 0.4475*x, which is at level 6
        // x is at level 5 (we switched the mod of x to level 5)
        evaluator.mod_switch_to_inplace(coeff1_plain, parms_ids[level+3]);
        evaluator.multiply_plain_inplace(x_encrypted, coeff1_plain);
        evaluator.rescale_to_next_inplace(x_encrypted);
        x_encrypted.scale() = scale;

        // To do the addition, we have to ensure all the terms have same parms_id and scale
        evaluator.mod_switch_to_inplace(x8_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x7_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x6_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x5_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x4_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x3_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[level+7]);
        evaluator.mod_switch_to_inplace(coeff0_plain, parms_ids[level+7]);

        // Calculate f(x)
        Ciphertext sign_encrypted;
        evaluator.add(x9_encrypted, x8_encrypted, sign_encrypted);
        evaluator.add_inplace(sign_encrypted, x7_encrypted);
        evaluator.add_inplace(sign_encrypted, x6_encrypted);
        evaluator.add_inplace(sign_encrypted, x5_encrypted);
        evaluator.add_inplace(sign_encrypted, x4_encrypted);
        evaluator.add_inplace(sign_encrypted, x3_encrypted);
        evaluator.add_inplace(sign_encrypted, x2_encrypted);
        evaluator.add_inplace(sign_encrypted, x_encrypted);
        evaluator.add_plain_inplace(sign_encrypted, coeff0_plain);

        // Calculate (m(t-1)-m(t))*(f(x)-1)
        Ciphertext result_encrypted;
        evaluator.mod_switch_to_inplace(recent_diff, parms_ids[level+7]);
        sign_encrypted.scale() = scale;
        recent_diff.scale() = scale;
        result_encrypted.scale() = scale;
        evaluator.multiply(sign_encrypted, recent_diff, result_encrypted);
        evaluator.relinearize_inplace(result_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(result_encrypted);
        result_encrypted.scale() = scale;

        decisions_encrypted.push_back(result_encrypted);
    }
    return decisions_encrypted;
}

// ***************************************************************************************************
// Data input/output functions
// ***************************************************************************************************

// Simulate the data import process. Return corresponding sample data according to the time
inline Ciphertext getSample(int time, int scale, CKKSEncoder &encoder, Encryptor &encryptor) {
    vector<double> prices = csv2vec("../data/apple_prices.csv");
    vector<double> sample;
    Plaintext sample_plain;
    Ciphertext sample_encrypted;
    sample.push_back(prices[time]);
    encoder.encode(sample, scale, sample_plain);
    encryptor.encrypt(sample_plain, sample_encrypted);
    return sample_encrypted;
}

// Append new sample data to the history data
inline void assembleSample(Ciphertext sample, vector<Ciphertext>& past_prices) {
    past_prices.push_back(sample);
}

// Decrytp a vector of Ciphertexts into a vector of doubles
inline vector<double> decryptVec(vector<Ciphertext>& ct, Decryptor &decryptor, CKKSEncoder &encoder) {
    vector<double> res;
    for (int i = 0; i < ct.size(); i++) {
        vector<double> result;
        Plaintext result_plain; 
        decryptor.decrypt(ct[i], result_plain);
        encoder.decode(result_plain, result);
        res.push_back(result[0]);
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
int main()
{
    // Set up the CKKS scheme.
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 16384;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 40, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40 }));
    double scale = pow(2.0, 30);
    int max_level = 12;

    // Set up Context
    auto context = SEALContext::Create(parms);
    print_parameters(context);
    cout << endl;

    // Generate public and private keys
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();
    auto relin_keys = keygen.relin_keys();
    GaloisKeys gal_keys = keygen.galois_keys();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;
    cout << endl;

    // Construct an array of parameters ids
    int p_i = 0;
    parms_id_type parms_ids[max_level];
    auto context_data = context->first_context_data();
    while (context_data)
    {
        parms_ids[p_i] = context_data->parms_id();
        context_data = context_data->next_context_data();
        p_i++;
    }

    // Sample data size
    int time_max = 200;

    cout << "Data Import Begins" << endl;
    auto import_start = high_resolution_clock::now();

    vector<Ciphertext> prices_encrypted;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample(i, scale, encoder, encryptor);
        assembleSample(sample, prices_encrypted);
    }

    auto import_stop = high_resolution_clock::now();
    auto import_duration = duration_cast<seconds>(import_stop - import_start); 
    cout << "Data Import Finished" << endl;
    cout << "Data Import Duration:  " << import_duration.count() << " seconds" << endl;
	cout << "\n";



	cout << "MACD Analysis Begins" << endl;
    auto macd_start = high_resolution_clock::now();

    // level 0: prices
    // level 1: wma12, wma26, wma_diff
    // level 2: wma9, macd
    vector<Ciphertext> wma12_encrypted = wma(prices_encrypted, 12, 0, encoder, encryptor, evaluator, scale, parms_ids);
    vector<Ciphertext> wma12_encrypted_sliced = slice(wma12_encrypted, 14, wma12_encrypted.size());

    vector<Ciphertext> wma26_encrypted = wma(prices_encrypted, 26, 0, encoder, encryptor, evaluator, scale, parms_ids);
    
    vector<Ciphertext> wma_diff_encrypted;
    for (int i = 0; i < wma26_encrypted.size(); i++) {
        Ciphertext tmp_diff;
        evaluator.sub(wma12_encrypted_sliced[i], wma26_encrypted[i], tmp_diff);
        wma_diff_encrypted.push_back(tmp_diff);
    }

    vector<Ciphertext> wma9_encrypted = wma(wma_diff_encrypted, 9, 1, encoder, encryptor, evaluator, scale, parms_ids);

    vector<Ciphertext> wma_diff_sliced = slice(wma_diff_encrypted, 9, wma_diff_encrypted.size());

    vector<Ciphertext> macd_encrypted;
    for (int i = 0; i < wma9_encrypted.size(); i++) {
        Ciphertext tmp_diff, tmp_macd;
        tmp_diff = wma_diff_sliced[i];
        evaluator.mod_switch_to_inplace(tmp_diff, parms_ids[2]);
        evaluator.sub(tmp_diff, wma9_encrypted[i], tmp_macd);
        macd_encrypted.push_back(tmp_macd);
    }

    auto macd_stop = high_resolution_clock::now();
    auto macd_duration = duration_cast<seconds>(macd_stop - macd_start); 
	cout << "MACD Analysis Finished" << endl;
    cout << "MACD Analysis Duration:  " << macd_duration.count() << " seconds" << endl;
	cout << "\n";

    cout << "Decision Analysis Begins" << endl;
    auto decision_start = high_resolution_clock::now();

    // level 2: m(t-1), m(t), m(t-1)-m(t)
    // level 3: norm*m(t), norm*m(t-1) // To normalize MACD signal
    // level 4: x = m(t-1)*m(t) // res_mult
    // level 5: x^2
    // level 6: x^3, x^4
    // level 7: x^5, x^6, x^7, x^8
    // level 8: x^9
    // level 9: f(x)
    // level 10: decision = (m(t)-m(t-1))*(f(x)-1)

    vector<Ciphertext> decisions_encrypted = decision(macd_encrypted, 2, 0.5, encoder, encryptor, 
        evaluator, relin_keys, scale, parms_ids);

    auto decision_stop = high_resolution_clock::now();
    auto decision_duration = duration_cast<seconds>(decision_stop - decision_start);
    cout << "Decision Analysis Finished" << endl;
    cout << "Decision Analysis Duration:  " << decision_duration.count() << " seconds" << endl;
    cout << "\n";

    cout << "Data Export Begins" << endl;
    auto export_start = high_resolution_clock::now();
    
    vector<double> wma12 = decryptVec(wma12_encrypted_sliced, decryptor, encoder);
    vector<double> wma26 = decryptVec(wma26_encrypted, decryptor, encoder);
    vector<double> wma_diff = decryptVec(wma_diff_sliced, decryptor, encoder);
    vector<double> wma9 = decryptVec(wma9_encrypted, decryptor, encoder);
    vector<double> macd = decryptVec(macd_encrypted, decryptor, encoder);
    vector<double> decisions = decryptVec(decisions_encrypted, decryptor, encoder);

    exportData(wma12, "../data/wma12_seal.csv");
    exportData(wma26, "../data/wma26_seal.csv");
    exportData(wma_diff, "../data/wma_diff_seal.csv");
    exportData(wma9, "../data/wma9_seal.csv");
    exportData(macd, "../data/macd_seal.csv");
    exportData(decisions, "../data/decisions_seal.csv");

    auto export_stop = high_resolution_clock::now();
    auto export_duration = duration_cast<seconds>(export_stop - export_start);
	cout << "Data Export Finished" << endl;
    cout << "Data Export Duration:  " << export_duration.count() << " seconds" << endl;
    cout << "\n";

    return 0;
}
