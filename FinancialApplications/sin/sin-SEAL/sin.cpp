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
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <stdlib.h>
#include <chrono>
#include <tuple>
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
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
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
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

// Return the sum of all the Ciphertexts in a vector
inline Ciphertext sum(vector<Ciphertext>& vec, Evaluator &evaluator, double scale) {
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


// Simulate the data import process. Return corresponding sample data according to the time
inline Ciphertext getSample(string file, int time, int scale, CKKSEncoder &encoder, Encryptor &encryptor) {
    vector<double> x = csv2vec(file);
    vector<double> sample;
    Plaintext sample_plain;
    Ciphertext sample_encrypted;
    sample.push_back(x[time]);
    encoder.encode(sample, scale, sample_plain);
    encryptor.encrypt(sample_plain, sample_encrypted);
    return sample_encrypted;
}


// Append new sample data to the history data
inline void assembleSample(Ciphertext sample, vector<Ciphertext>& samples) {
    samples.push_back(sample);
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


// Calculate the percentage error between two lists of data
inline double pe(vector<double> ls1, vector<double> ls2) {
    if (ls1.size() != ls2.size()) {
        cout << "Error: the input lists are of different length!" << endl;
        return 0.0;
    } else {
        int n = ls1.size();
        double sm = 0.0;
        for (int i = 0; i < n; i++) {
            sm += abs((ls1[i]-ls2[i])/ls2[i])*100;
        }
        return sm/n;
    }
}


// Get the i-th column of a matrix
inline vector<double> get_column(vector<vector<double>> matrix, int index) {
    vector<double> col;
    for (int i = 0; i < matrix.size(); i++) {
        col.push_back(matrix[i][index]);
    }
    return col;
}


// Calculate the factorial of a number
inline long factorial(int x) {
    long res = 1;
    for (int i = 1; i <= x; ++i) {
        res *= i;
    }
    return res;
}


inline Ciphertext s0_sin(Ciphertext x, int r, int terms, CKKSEncoder &encoder, Evaluator &evaluator, 
        Decryptor &decryptor, RelinKeys &relin_keys, double scale) {
    Ciphertext result, x_encrypted_coeff0, x2_encrypted, x3_encrypted, x3_encrypted_coeff1, x5_encrypted, x5_encrypted_coeff2;

    // Encrypt parameters
    Plaintext coeff0, coeff1, coeff2;
    encoder.encode(0.125, scale, coeff0);
    encoder.encode(-0.000325520833, scale, coeff1);
    encoder.encode(0.0000002543131504167, scale, coeff2);

    // Calculate 0.125*x
    evaluator.multiply_plain(x, coeff0, x_encrypted_coeff0);
    evaluator.rescale_to_next_inplace(x_encrypted_coeff0);
    x_encrypted_coeff0.scale() = scale;

    // Calculate x^2
    evaluator.square(x, x2_encrypted);
    evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x2_encrypted);
    x2_encrypted.scale() = scale;

    // Calculate -0.000325520833*x^3
    evaluator.mod_switch_to_inplace(x, x2_encrypted.parms_id());
    evaluator.multiply(x, x2_encrypted, x3_encrypted);
    evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x3_encrypted);
    x3_encrypted.scale() = scale;
    evaluator.mod_switch_to_inplace(coeff1, x3_encrypted.parms_id());
    evaluator.multiply_plain(x3_encrypted, coeff1, x3_encrypted_coeff1);
    evaluator.rescale_to_next_inplace(x3_encrypted_coeff1);
    x3_encrypted_coeff1.scale() = scale;

    // Calculate 0.0000002543131504167*x^5
    evaluator.mod_switch_to_inplace(x2_encrypted, x3_encrypted.parms_id());
    evaluator.multiply(x2_encrypted, x3_encrypted, x5_encrypted);
    evaluator.relinearize_inplace(x5_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x5_encrypted);
    x5_encrypted.scale() = scale;
    evaluator.mod_switch_to_inplace(coeff2, x5_encrypted.parms_id());
    evaluator.multiply_plain(x5_encrypted, coeff2, x5_encrypted_coeff2);
    evaluator.rescale_to_next_inplace(x5_encrypted_coeff2);
    x5_encrypted_coeff2.scale() = scale;

    // Calculate s0(x)
    evaluator.mod_switch_to_inplace(x_encrypted_coeff0, x3_encrypted_coeff1.parms_id());
    evaluator.add_inplace(x3_encrypted_coeff1, x_encrypted_coeff0);
    evaluator.mod_switch_to_inplace(x3_encrypted_coeff1, x5_encrypted_coeff2.parms_id());
    evaluator.add(x5_encrypted_coeff2, x3_encrypted_coeff1, result);

    return result;
}


inline Ciphertext c0(Ciphertext x, int r, int terms, CKKSEncoder &encoder, Evaluator &evaluator, 
        Decryptor &decryptor, RelinKeys &relin_keys, double scale) {
    Ciphertext result, x2_encrypted, x2_encrypted_coeff1, x4_encrypted, x4_encrypted_coeff2;

    // Encrypt parameters
    Plaintext coeff0, coeff1, coeff2;
    encoder.encode(0.125, scale, coeff0);
    encoder.encode(-0.000325520833, scale, coeff1);
    encoder.encode(0.0000002543131504167, scale, coeff2);

    // Calculate -0.0078125*x^2
    evaluator.square(x, x2_encrypted);
    evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x2_encrypted);
    x2_encrypted.scale() = scale;
    evaluator.mod_switch_to_inplace(coeff1, x2_encrypted.parms_id());
    evaluator.multiply_plain(x2_encrypted, coeff1, x2_encrypted_coeff1);
    evaluator.rescale_to_next_inplace(x2_encrypted_coeff1);
    x2_encrypted_coeff1.scale() = scale;

    // Calculate 1.017252604167*x^4
    evaluator.square(x2_encrypted, x4_encrypted);
    evaluator.relinearize_inplace(x4_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x4_encrypted);
    x4_encrypted.scale() = scale;
    evaluator.mod_switch_to_inplace(coeff2, x4_encrypted.parms_id());
    evaluator.multiply_plain(x4_encrypted, coeff2, x4_encrypted_coeff2);
    evaluator.rescale_to_next_inplace(x4_encrypted_coeff2);
    x4_encrypted_coeff2.scale() = scale;

    // Calculate c0(x)
    evaluator.mod_switch_to_inplace(x2_encrypted_coeff1, x4_encrypted_coeff2.parms_id());
    evaluator.add_inplace(x4_encrypted_coeff2, x2_encrypted_coeff1);
    evaluator.mod_switch_to_inplace(coeff0, x4_encrypted_coeff2.parms_id());
    evaluator.add_plain(x4_encrypted_coeff2, coeff0, result);

    return result;
}


int main() {
    // SEAL settings
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 32768;

    cout << "max coeff_modulus bit-length: " << CoeffModulus::MaxBitCount(poly_modulus_degree) << endl;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {40, 30, 30, 30, 30, 30, 
    //         30, 30, 30, 30, 30, 40}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 30, 30, 30, 
            30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50}));
    double scale = pow(2.0, 30);
    SEALContext context(parms);
    print_parameters(context);
    cout << endl;

    // Generate public and private keys
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;
    cout << endl;

    int r = 3;
    int terms = 3;
    int time_max = 158;

    cout << "Data Import Starts" << endl;
    vector<Ciphertext> x_samples;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample("../../data/pi_samples.csv", i, scale, encoder, encryptor);
        assembleSample(sample, x_samples);
    }
    cout << "Data Import Ends" << endl;

    // Apply recursive PID on Bergman Minimal Model
    // Only the last two errors and the previous sum are passed into the algorithm
    cout << "Sine Interpolation Test Begins" << endl;
    auto sin_start = high_resolution_clock::now();

    vector<Ciphertext> sin_encrypted;
    
    for (int i = 0; i < x_samples.size(); i++) {

        Ciphertext sk0, ck0, x;
        Plaintext coeff;
        x = x_samples[i];

        // cout << "Start Calculating s0 and c0" << endl;

        encoder.encode((2.0), scale, coeff);
        sk0 = s0_sin(x, r, terms, encoder, evaluator, decryptor, relin_keys, scale);
        ck0 = c0(x, r, terms, encoder, evaluator, decryptor, relin_keys, scale);

        // cout << "Finish Calculating s0 and c0" << endl;

        // Decrypt
        // vector<double> sk0_val;
        // Plaintext sk0_val_plain;
        // decryptor.decrypt(sk0, sk0_val_plain);
        // encoder.decode(sk0_val_plain, sk0_val);
        // cout << "sk0: " << sk0_val[0] << endl;


        // cout << "Start 1st Iteration" << endl;
        Ciphertext sk1, ck1, sk_square, ck_square;
        // Calculate ck1
        evaluator.square(sk0, sk_square);
        evaluator.relinearize_inplace(sk_square, relin_keys);
        evaluator.rescale_to_next_inplace(sk_square);
        sk_square.scale() = scale;
        evaluator.square(ck0, ck_square);
        evaluator.relinearize_inplace(ck_square, relin_keys);
        evaluator.rescale_to_next_inplace(ck_square);
        ck_square.scale() = scale;
        evaluator.mod_switch_to_inplace(ck_square, sk_square.parms_id());
        evaluator.sub(ck_square, sk_square, ck1);

        // Calculate sk1
        evaluator.mod_switch_to_inplace(ck0, sk0.parms_id());
        evaluator.multiply(sk0, ck0, sk1);
        evaluator.relinearize_inplace(sk1, relin_keys);
        evaluator.rescale_to_next_inplace(sk1);
        sk1.scale() = scale;
        evaluator.mod_switch_to_inplace(coeff, sk1.parms_id());
        evaluator.multiply_plain_inplace(sk1, coeff);
        evaluator.rescale_to_next_inplace(sk1);
        sk1.scale() = scale;

        // cout << "Finish 1st Iteration" << endl;

        // Decrypt
        // vector<double> sk1_val;
        // Plaintext sk1_val_plain;
        // decryptor.decrypt(sk1, sk1_val_plain);
        // encoder.decode(sk1_val_plain, sk1_val);
        // cout << "sk1: " << sk1_val[0] << endl;


        // cout << "Start 2nd Iteration" << endl;
        Ciphertext sk2, ck2, sk_square2, ck_square2;
        // Calculate ck2
        evaluator.square(sk1, sk_square2);
        evaluator.relinearize_inplace(sk_square2, relin_keys);
        evaluator.rescale_to_next_inplace(sk_square2);
        sk_square2.scale() = scale;
        evaluator.square(ck1, ck_square2);
        evaluator.relinearize_inplace(ck_square2, relin_keys);
        evaluator.rescale_to_next_inplace(ck_square2);
        ck_square2.scale() = scale;
        evaluator.mod_switch_to_inplace(ck_square2, sk_square2.parms_id());
        evaluator.sub(ck_square2, sk_square2, ck2);

        // Calculate sk2
        evaluator.mod_switch_to_inplace(ck1, sk1.parms_id());
        evaluator.multiply(sk1, ck1, sk2);
        evaluator.relinearize_inplace(sk2, relin_keys);
        evaluator.rescale_to_next_inplace(sk2);
        sk2.scale() = scale;
        evaluator.mod_switch_to_inplace(coeff, sk2.parms_id());
        evaluator.multiply_plain_inplace(sk2, coeff);
        evaluator.rescale_to_next_inplace(sk2);
        sk2.scale() = scale;

        // cout << "Finish 2nd Iteration" << endl;

        // Decrypt
        // vector<double> sk2_val;
        // Plaintext sk2_val_plain;
        // decryptor.decrypt(sk2, sk2_val_plain);
        // encoder.decode(sk2_val_plain, sk2_val);
        // cout << "sk2: " << sk2_val[0] << endl;


        // cout << "Start 3rd Iteration" << endl;
        Ciphertext sk3, ck3, sk_square3, ck_square3;
        // Calculate ck3
        evaluator.square(sk2, sk_square3);
        evaluator.relinearize_inplace(sk_square3, relin_keys);
        evaluator.rescale_to_next_inplace(sk_square3);
        sk_square3.scale() = scale;
        evaluator.square(ck2, ck_square3);
        evaluator.relinearize_inplace(ck_square3, relin_keys);
        evaluator.rescale_to_next_inplace(ck_square3);
        ck_square3.scale() = scale;
        evaluator.mod_switch_to_inplace(ck_square3, sk_square3.parms_id());
        evaluator.sub(ck_square3, sk_square3, ck3);

        // Calculate sk3
        evaluator.mod_switch_to_inplace(ck2, sk2.parms_id());
        evaluator.multiply(sk2, ck2, sk3);
        evaluator.relinearize_inplace(sk3, relin_keys);
        evaluator.rescale_to_next_inplace(sk3);
        sk3.scale() = scale;
        evaluator.mod_switch_to_inplace(coeff, sk3.parms_id());
        evaluator.multiply_plain_inplace(sk3, coeff);
        evaluator.rescale_to_next_inplace(sk3);
        sk3.scale() = scale;

        // cout << "Finish 3rd Iteration" << endl;

        // Decrypt
        // vector<double> sk3_val;
        // Plaintext sk3_val_plain;
        // decryptor.decrypt(sk3, sk3_val_plain);
        // encoder.decode(sk3_val_plain, sk3_val);
        // cout << "sk3: " << sk3_val[0] << endl;

        sin_encrypted.push_back(sk2);
    }

    cout << "Sine Interpolation Test Ends" << endl;
    auto sin_stop = high_resolution_clock::now();
    auto sin_duration = duration_cast<seconds>(sin_stop - sin_start);
    cout << "Sine Interpolation Test Duration: " << sin_duration.count() << " seconds" << endl;
    cout << endl;

    cout << "Output data" << endl;
    cout << endl;

    vector<double> sin_plain = decryptVec(sin_encrypted, decryptor, encoder);
    exportData(sin_plain, "../../data/sine_interpolated_seal_2.csv");
    cout << endl;
}
