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
#include <tuple>
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

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

inline Ciphertext pid_controller(vector<Ciphertext>& error_hist, CKKSEncoder &encoder,
        Evaluator &evaluator, Decryptor &decryptor, RelinKeys &relin_keys, double scale) {
        
    Plaintext u_steady, Kp, Ki_tor_i, Kd, neg_one;
    encoder.encode((16.67), scale, u_steady);
    encoder.encode((0.5), scale, Kp);
    encoder.encode((0.01), scale, Ki_tor_i);
    encoder.encode((0.05), scale, Kd);
    encoder.encode((-1), scale, neg_one);

    Ciphertext result_encrypted, last_signal, sec_last_signal, first_term, second_term, third_term, cur_sum, dif;

    last_signal = error_hist[error_hist.size()-1];
    sec_last_signal = error_hist[error_hist.size()-2];

    // Kp * error_hist[-1]
    evaluator.multiply_plain(last_signal, Kp, first_term);
    evaluator.rescale_to_next_inplace(first_term);
    first_term.scale() = scale;

    // Ki * sum(error_hist) / tor_i
    cur_sum = sum(error_hist, evaluator, scale);
    evaluator.multiply_plain(cur_sum, Ki_tor_i, second_term);
    evaluator.rescale_to_next_inplace(second_term);
    second_term.scale() = scale;

    // Kd * (error_hist[-1] - error_hist[-2])
    evaluator.multiply_plain(sec_last_signal, neg_one, third_term);
    evaluator.rescale_to_next_inplace(third_term);
    third_term.scale() = scale;
    evaluator.mod_switch_to_inplace(last_signal, third_term.parms_id());
    evaluator.add_inplace(third_term, last_signal);
    evaluator.mod_switch_to_inplace(Kd, third_term.parms_id());
    evaluator.multiply_plain_inplace(third_term, Kd);
    evaluator.rescale_to_next_inplace(third_term);
    third_term.scale() = scale;

    // Add all
    evaluator.add(first_term, second_term, result_encrypted);
    evaluator.mod_switch_to_inplace(result_encrypted, third_term.parms_id());
    evaluator.add_inplace(result_encrypted, third_term);
    evaluator.mod_switch_to_inplace(u_steady, result_encrypted.parms_id());
    evaluator.add_plain_inplace(result_encrypted, u_steady);

    return result_encrypted;
}

inline tuple<Ciphertext, Ciphertext> pid_controller_recursive(vector<Ciphertext>& error_hist, CKKSEncoder &encoder,
        Evaluator &evaluator, Decryptor &decryptor, RelinKeys &relin_keys, double scale, Ciphertext prev_sum) {
        
    Plaintext u_steady, Kp, Ki_tor_i, Kd, neg_one;
    encoder.encode((16.67), scale, u_steady);
    encoder.encode((0.5), scale, Kp);
    encoder.encode((0.01), scale, Ki_tor_i);
    encoder.encode((0.05), scale, Kd);
    encoder.encode((-1), scale, neg_one);

    Ciphertext result_encrypted, last_signal, sec_last_signal, first_term, second_term, third_term, cur_sum, dif;

    last_signal = error_hist[error_hist.size()-1];
    sec_last_signal = error_hist[error_hist.size()-2];

    // Kp * error_hist[-1]
    evaluator.multiply_plain(last_signal, Kp, first_term);
    evaluator.rescale_to_next_inplace(first_term);
    first_term.scale() = scale;

    // Ki * sum(error_hist) / tor_i
    evaluator.mod_switch_to_inplace(prev_sum, last_signal.parms_id());
    evaluator.add(last_signal, prev_sum, cur_sum);
    evaluator.multiply_plain(cur_sum, Ki_tor_i, second_term);
    evaluator.rescale_to_next_inplace(second_term);
    second_term.scale() = scale;

    // Kd * (error_hist[-1] - error_hist[-2])
    evaluator.multiply_plain(sec_last_signal, neg_one, third_term);
    evaluator.rescale_to_next_inplace(third_term);
    third_term.scale() = scale;
    evaluator.mod_switch_to_inplace(last_signal, third_term.parms_id());
    evaluator.add_inplace(third_term, last_signal);
    evaluator.mod_switch_to_inplace(Kd, third_term.parms_id());
    evaluator.multiply_plain_inplace(third_term, Kd);
    evaluator.rescale_to_next_inplace(third_term);
    third_term.scale() = scale;

    // Add all
    evaluator.add(first_term, second_term, result_encrypted);
    evaluator.mod_switch_to_inplace(result_encrypted, third_term.parms_id());
    evaluator.add_inplace(result_encrypted, third_term);
    evaluator.mod_switch_to_inplace(u_steady, result_encrypted.parms_id());
    evaluator.add_plain_inplace(result_encrypted, u_steady);

    return make_tuple(result_encrypted, cur_sum);
}

int main() {
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    double scale = pow(2.0, 40);
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

    cout << "Reading input data" << endl;

    vector<Ciphertext> error_hist;
    vector<double> input = csv2vec("../data/white_noise.csv");

    // Encrypt input
    for (int i = 0; i < input.size(); i++) {
        Plaintext tmp_plain;
        encoder.encode(input[i], scale, tmp_plain);
        Ciphertext tmp_encrypted;
        encryptor.encrypt(tmp_plain, tmp_encrypted);
        error_hist.push_back(tmp_encrypted);
    }

    cout << "Applying PID controller" << endl;

    // PID
    vector<Ciphertext> result_encrypted;
    for (int i = 2; i < error_hist.size(); i++) {
        vector<Ciphertext> sliced_hist = slice(error_hist, 0, i);
        Ciphertext tmp = pid_controller(sliced_hist, encoder, evaluator, decryptor, relin_keys, scale);
        result_encrypted.push_back(tmp);
        
    }

    cout << "Applying recursive PID controller" << endl;

    vector<Ciphertext> result_rec_encrypted;
    vector<Ciphertext> prev_sum_vec;
    // init
    Plaintext zero;
    encoder.encode((0), scale, zero);
    Ciphertext zero_encrypted;
    encryptor.encrypt(zero, zero_encrypted);
    vector<Ciphertext> sliced_hist = slice(error_hist, 0, 2);
    tuple<Ciphertext, Ciphertext> tup = pid_controller_recursive(sliced_hist, encoder, evaluator, decryptor, relin_keys, scale, zero_encrypted);
    result_rec_encrypted.push_back(get<0>(tup));
    prev_sum_vec.push_back(get<1>(tup));

    for (int i = 3; i < error_hist.size(); i++) {
        vector<Ciphertext> sliced_hist = slice(error_hist, i-2, i);
        tuple<Ciphertext, Ciphertext> tup = pid_controller_recursive(sliced_hist, encoder, evaluator, decryptor, 
            relin_keys, scale, prev_sum_vec[prev_sum_vec.size()-1]);
        result_rec_encrypted.push_back(get<0>(tup));
        prev_sum_vec.push_back(get<1>(tup));
    }

    cout << "Output data" << endl;

    vector<double> result_plaintext = decryptVec(result_encrypted, decryptor, encoder);
    exportData(result_plaintext, "../data/seal_result.csv");

    vector<double> result_rec_plaintext = decryptVec(result_rec_encrypted, decryptor, encoder);
    exportData(result_rec_plaintext, "../data/seal_result_recursive.csv");
}