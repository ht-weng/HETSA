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
#include "seal/seal.h"

using namespace std;
using namespace seal;

//********************************************************************************
// SEAL helper functions
//********************************************************************************
/*
Helper function: Prints the name of the example in a fancy banner.
*/
inline void print_example_banner(std::string title)
{
    if (!title.empty())
    {
        std::size_t title_length = title.length();
        std::size_t banner_length = title_length + 2 * 10;
        std::string banner_top = "+" + std::string(banner_length - 2, '-') + "+";
        std::string banner_middle =
            "|" + std::string(9, ' ') + title + std::string(9, ' ') + "|";

        std::cout << std::endl
            << banner_top << std::endl
            << banner_middle << std::endl
            << banner_top << std::endl;
    }
}

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
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    stream << std::hex << std::setfill('0')
        << std::setw(16) << parms_id[0] << " "
        << std::setw(16) << parms_id[1] << " "
        << std::setw(16) << parms_id[2] << " "
        << std::setw(16) << parms_id[3] << " ";

    /*
    Restore the old std::cout formatting.
    */
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

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}


/*
Helper function: Prints a matrix of values.
*/
template<typename T>
inline void print_matrix(std::vector<T> matrix, std::size_t row_size)
{
    /*
    We're not going to print every column of the matrix (there are 2048). Instead
    print this many slots from beginning and end of the matrix.
    */
    std::size_t print_size = 5;

    std::cout << std::endl;
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++)
    {
        std::cout << std::setw(3) << std::right << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = row_size - print_size; i < row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i]
            << ((i != row_size - 1) ? "," : " ]\n");
    }
    std::cout << "    [";
    for (std::size_t i = row_size; i < row_size + print_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i]
            << ((i != 2 * row_size - 1) ? "," : " ]\n");
    }
    std::cout << std::endl;
};

/*
Helper function: Print line number.
*/
inline void print_line(int line_number)
{
    std::cout << "Line " << std::setw(3) << line_number << " --> ";
}

//////////////////////////////////////////////////////////////////////////////////
// MACD helper functions
//////////////////////////////////////////////////////////////////////////////////
inline Ciphertext sum(vector<Ciphertext>& v, CKKSEncoder &encoder, Encryptor &encryptor, 
    Evaluator &evaluator, double scale) {
    Ciphertext result_encrypted = v[0];

    for(int i = 1; i < v.size(); i++) {
        Ciphertext val = v[i];
        evaluator.add_inplace(result_encrypted, val);
    }
    return result_encrypted;
}

inline vector<Ciphertext> slice(vector<Ciphertext>& v, int start=0, int end=-1) {
    int oldlen = v.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<Ciphertext> nv;

    for (int i = 0; i < newlen; i++) {
        nv.push_back(v[start+i]);
    }
    return nv;
}

inline void print_vector(const vector<double>& v) {
    for(int i=0; i < v.size(); i++) {
        cout << v[i] << ' ';
    }
    cout << "\n";
}

inline vector<double> csv2vec(string inputFileName) {
    vector<double> data;
    ifstream inputFile(inputFileName);
    int l = 0;

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
            cout << "NaN found in file " << inputFileName << " line " << l
                 << endl;
            e.what();
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}

inline vector<Ciphertext> wma(vector<Ciphertext>& s, int n, int level, CKKSEncoder &encoder, 
    Encryptor &encryptor, Evaluator &evaluator, double scale, parms_id_type *parms_ids) {
    vector<Ciphertext> wma;
    vector<Plaintext> weights;

    // generate a list of weights of the window size
    for (int i = 0; i < n; i++) {
        vector<double> w;
        Plaintext w_plain;
        w.push_back(2.0*(i+1.0)/(n*(n+1.0)));
        encoder.encode(w, scale, w_plain);
        weights.push_back(w_plain);
    }
    // multiply corresponding data points and weights to get WMA
    for (int i = 0; i < s.size()-n; i++) {
        vector<Ciphertext> s_sliced;
        vector<Ciphertext> window;
        s_sliced = slice(s, i, i+n);

        for (int j = 0; j < n; j++) {
            Ciphertext tmp = s_sliced[j];
            Plaintext tmp_weight = weights[j];
            evaluator.mod_switch_to_inplace(tmp_weight, parms_ids[level]);
            evaluator.multiply_plain_inplace(tmp, tmp_weight);
            evaluator.rescale_to_next_inplace(tmp);
            tmp.scale() = scale;
            window.push_back(tmp);
        }

        Ciphertext sumup = sum(window, encoder, encryptor, evaluator, scale);
        wma.push_back(sumup);
    }
    return wma;
}

inline Ciphertext getSample(int time, int scale, CKKSEncoder &encoder, Encryptor &encryptor) {
    vector<double> prices = csv2vec("apple_prices.csv");
    vector<double> sample;
    Plaintext sample_plain;
    Ciphertext sample_encrypted;
    sample.push_back(prices[time]);
    encoder.encode(sample, scale, sample_plain);
    encryptor.encrypt(sample_plain, sample_encrypted);
    return sample_encrypted;
}

inline void assembleSample(Ciphertext sample, vector<Ciphertext>& past_prices) {
    past_prices.push_back(sample);
}

//********************************************************************************
// Main function
//********************************************************************************
int main()
{
    //////////////////////////////////////////////////////////////////////////////
    // SEAL settings
    //////////////////////////////////////////////////////////////////////////////
    // Set up the CKKS scheme.
    EncryptionParameters parms(scheme_type::CKKS);
    // macd and tanh
    // size_t poly_modulus_degree = 8192;
    size_t poly_modulus_degree = 16384;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    // macd setting
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 40, 30, 30, 30, 40 }));
    // tanh setting
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 30, 25, 25, 25, 25, 25, 25, 25, 30 }));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 40, 30, 30, 30, 30, 30, 30, 30, 30, 30, 40 }));
    double scale = pow(2.0, 30);
    int max_level = 11;

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

    //////////////////////////////////////////////////////////////////////////////
    // Simulation of Receiving Data From Client
    //////////////////////////////////////////////////////////////////////////////
    int time_max = 500;

    vector<Ciphertext> prices_encrypted;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample(i, scale, encoder, encryptor);
        assembleSample(sample, prices_encrypted);
    }

    //////////////////////////////////////////////////////////////////////////////
    // MACD
    //////////////////////////////////////////////////////////////////////////////
    // level 0: prices
    // level 1: wma12, wma26, wma_diff
    // level 2: wma9, macd

    vector<Ciphertext> wma12_encrypted;
    wma12_encrypted = wma(prices_encrypted, 12, 0, encoder, encryptor, evaluator, scale, parms_ids);
    vector<Ciphertext> wma12_encrypted_sliced = slice(wma12_encrypted, 14, wma12_encrypted.size());

    vector<Ciphertext> wma26_encrypted;
    wma26_encrypted = wma(prices_encrypted, 26, 0, encoder, encryptor, evaluator, scale, parms_ids);
    
    vector<Ciphertext> wma_diff_encrypted;
    for (int i = 0; i < wma26_encrypted.size(); i++) {
        Ciphertext tmp_diff;
        evaluator.sub(wma12_encrypted_sliced[i], wma26_encrypted[i], tmp_diff);
        wma_diff_encrypted.push_back(tmp_diff);
    }

    vector<Ciphertext> wma9_encrypted;
    wma9_encrypted = wma(wma_diff_encrypted, 9, 1, encoder, encryptor, evaluator, scale, parms_ids);

    vector<Ciphertext> wma_diff_sliced = slice(wma_diff_encrypted, 9, wma_diff_encrypted.size());

    vector<Ciphertext> macd_encrypted;
    for (int i = 0; i < wma9_encrypted.size(); i++) {
        Ciphertext tmp_diff, tmp_macd;
        tmp_diff = wma_diff_sliced[i];
        evaluator.mod_switch_to_inplace(tmp_diff, parms_ids[2]);
        evaluator.sub(tmp_diff, wma9_encrypted[i], tmp_macd);
        macd_encrypted.push_back(tmp_macd);
    }

    cout << "MACD signals calculation done" << endl;

    //////////////////////////////////////////////////////////////////////////////
    // Decision
    //////////////////////////////////////////////////////////////////////////////
    // level 2: m(t), m(t-1), m(t)-m(t-1)
    // level 3: 0.25*m(t), 0.25*m(t-1) // To normalize MACD signal
    // level 4: x = m(t)*m(t-1) // result_mult
    // level 5: x^2
    // level 6: x^3 = x^2*x
    // level 7: x^5 = x^3*x^2
    // level 8: f(x) = 0.375*x^5 - 1.25*x^3 + 1.875*x
    // level 9: decision = (m(t)-m(t-1))*f(x)

    vector<Ciphertext> decisions_encrypted;
    for (int i = 1; i < macd_encrypted.size(); i++) {

        // First store the latest two macd signals m(t) and m(t-1)
        Ciphertext mt = macd_encrypted[i];
        Ciphertext mt_1 = macd_encrypted[i-1];

        Plaintext coeff_norm_plain;
        encoder.encode(0.25, scale, coeff_norm_plain);
        evaluator.mod_switch_to_inplace(coeff_norm_plain, parms_ids[2]);

        evaluator.multiply_plain_inplace(mt, coeff_norm_plain);
        evaluator.rescale_to_next_inplace(mt);
        mt.scale() = scale;
        evaluator.multiply_plain_inplace(mt_1, coeff_norm_plain);
        evaluator.rescale_to_next_inplace(mt_1);
        mt_1.scale() = scale;

        Ciphertext recent_diff;
        evaluator.sub(mt, mt_1, recent_diff);
        Ciphertext recent_mult;
        // recent_mult is at level 4
        evaluator.multiply(mt, mt_1, recent_mult);
        evaluator.relinearize_inplace(recent_mult, relin_keys);
        evaluator.rescale_to_next_inplace(recent_mult);
        recent_mult.scale() = scale;

        // Set up coefficients
        Plaintext coeff1_plain, coeff2_plain, coeff3_plain, offset_plain;
        Ciphertext offset_encrypted;
        encoder.encode(0.375, scale, coeff1_plain);
        encoder.encode((-1.25), scale, coeff2_plain);
        encoder.encode(1.875, scale, coeff3_plain);
        encoder.encode(1.0, scale, offset_plain);
        encryptor.encrypt(offset_plain, offset_encrypted);

        // Calculate x^2, which is at level 5
        Ciphertext x2_encrypted;
        evaluator.square(recent_mult, x2_encrypted);
        evaluator.relinearize_inplace(x2_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x2_encrypted);
        x2_encrypted.scale() = scale;

        // Calculate x^3, which is at level 6
        Ciphertext x3_encrypted;
        // Switch mod to ensure parameters id match
        // x^2 is at level 5
        evaluator.mod_switch_to_inplace(recent_mult, parms_ids[5]);
        evaluator.multiply(recent_mult, x2_encrypted, x3_encrypted);
        evaluator.relinearize_inplace(x3_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x3_encrypted);
        x3_encrypted.scale() = scale;

        // Calculate x^5, which is at level 7
        Ciphertext x5_encrypted;
        // Switch mod to ensure parameters id match
        // x^3 is at level 6
        evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[6]);
        evaluator.multiply(x2_encrypted, x3_encrypted, x5_encrypted);
        evaluator.relinearize_inplace(x5_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(x5_encrypted);
        x5_encrypted.scale() = scale;

        // Calculate 0.375*x^5, which is at level 8
        // Switch mod to ensure parameters id match
        // x^5 is at level 7
        evaluator.mod_switch_to_inplace(coeff1_plain, parms_ids[7]);
        evaluator.multiply_plain_inplace(x5_encrypted, coeff1_plain);
        evaluator.rescale_to_next_inplace(x5_encrypted);
        x5_encrypted.scale() = scale;

        // Calculate -1.25*x^3, which is at level 7
        // Switch mod to ensure parameters id match
        // x^3 is at level 6
        evaluator.mod_switch_to_inplace(coeff2_plain, parms_ids[6]);
        evaluator.multiply_plain_inplace(x3_encrypted, coeff2_plain);
        evaluator.rescale_to_next_inplace(x3_encrypted);
        x3_encrypted.scale() = scale;

        // Calculate 1.875*x
        // Switch mod to ensure parameters id match
        // Note that previously we haved switched recent_mult' level to 5
        evaluator.mod_switch_to_inplace(coeff3_plain, parms_ids[5]);
        evaluator.multiply_plain_inplace(recent_mult, coeff3_plain);
        evaluator.rescale_to_next_inplace(recent_mult);
        recent_mult.scale() = scale;

        // To do the addition, we have to ensure that the terms have same parms_id and scale
        evaluator.mod_switch_to_inplace(x3_encrypted, parms_ids[8]);
        evaluator.mod_switch_to_inplace(recent_mult, parms_ids[8]);
        evaluator.mod_switch_to_inplace(offset_encrypted, parms_ids[8]);

        Ciphertext sign_encrypted;
        evaluator.add(x5_encrypted, x3_encrypted, sign_encrypted);
        evaluator.add_inplace(sign_encrypted, recent_mult);
        evaluator.add_inplace(sign_encrypted, offset_encrypted);

        Ciphertext result_encrypted;
        evaluator.mod_switch_to_inplace(recent_diff, parms_ids[8]);
        sign_encrypted.scale() = scale;
        recent_diff.scale() = scale;
        result_encrypted.scale() = scale;
        evaluator.multiply(sign_encrypted, recent_diff, result_encrypted);
        evaluator.relinearize_inplace(result_encrypted, relin_keys);
        evaluator.rescale_to_next_inplace(result_encrypted);
        result_encrypted.scale() = scale;

        decisions_encrypted.push_back(result_encrypted);
    }

    vector<double> decisions;
    for (int i = 0; i < decisions_encrypted.size(); i++) {
        vector<double> result;
        Plaintext result_plain; 
        decryptor.decrypt(decisions_encrypted[i], result_plain);
        encoder.decode(result_plain, result);
        decisions.push_back(result[0]);
    }

    cout << "Decisions calculation done" << endl;
    // //////////////////////////////////////////////////////////////////////////////
    // // tanh
    // //////////////////////////////////////////////////////////////////////////////
    // // macd_encrypted is at level 3
    // vector<Ciphertext> tanh_encrypted;
    // for (int i = 0; i < macd_encrypted.size(); i++) {
    //     // Set up coefficients
    //     Plaintext coeff1_plain, coeff2_plain, coeff3_plain;
    //     encoder.encode(0.375, scale, coeff1_plain);
    //     encoder.encode((-1.25), scale, coeff2_plain);
    //     encoder.encode(1.875, scale, coeff3_plain);

    //     Ciphertext x_encrypted = macd_encrypted[i];

    //     // Calculate x^2, which is at level 1
    //     Ciphertext x2_encrypted;
    //     evaluator.square(x_encrypted, x2_encrypted);
    //     evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    //     evaluator.rescale_to_next_inplace(x2_encrypted);
    //     x2_encrypted.scale() = scale;

    //     // Calculate x^3, which is at level 2
    //     Ciphertext x3_encrypted;
    //     // Switch mod to ensure parameters id match
    //     evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[3]);
    //     evaluator.multiply(x_encrypted, x2_encrypted, x3_encrypted);
    //     evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    //     evaluator.rescale_to_next_inplace(x3_encrypted);
    //     x3_encrypted.scale() = scale;

    //     // Calculate x^5, which is at level 3
    //     Ciphertext x5_encrypted;
    //     // Switch mod to ensure parameters id match
    //     evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[4]);
    //     evaluator.multiply(x2_encrypted, x3_encrypted, x5_encrypted);
    //     evaluator.relinearize_inplace(x5_encrypted, relin_keys);
    //     evaluator.rescale_to_next_inplace(x5_encrypted);
    //     x5_encrypted.scale() = scale;

    //     // Calculate 0.375*x^5, which is at level 4
    //     // Switch mod to ensure parameters id match
    //     evaluator.mod_switch_to_inplace(coeff1_plain, parms_ids[5]);
    //     evaluator.multiply_plain_inplace(x5_encrypted, coeff1_plain);
    //     evaluator.rescale_to_next_inplace(x5_encrypted);
    //     x5_encrypted.scale() = scale;

    //     // Calculate -1.25*x^3, which is at level 3
    //     // Switch mod to ensure parameters id match
    //     evaluator.mod_switch_to_inplace(coeff2_plain, parms_ids[4]);
    //     evaluator.multiply_plain_inplace(x3_encrypted, coeff2_plain);
    //     evaluator.rescale_to_next_inplace(x3_encrypted);
    //     x3_encrypted.scale() = scale;

    //     // Calculate 1.875*x, which is at level 2
    //     // Switch mod to ensure parameters id match
    //     // Note that previously we haved switched x_encrypted' level to 1
    //     evaluator.mod_switch_to_inplace(coeff3_plain, parms_ids[3]);
    //     evaluator.multiply_plain_inplace(x_encrypted, coeff3_plain);
    //     evaluator.rescale_to_next_inplace(x_encrypted);
    //     x_encrypted.scale() = scale;

    //     // To do the addition, we have to ensure that the terms have same parms_id and scale
    //     evaluator.mod_switch_to_inplace(x3_encrypted, parms_ids[6]);
    //     evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[6]);

    //     Ciphertext result_encrypted;
    //     evaluator.add(x5_encrypted, x3_encrypted, result_encrypted);
    //     evaluator.add_inplace(result_encrypted, x_encrypted);

    //     tanh_encrypted.push_back(result_encrypted);
    // }

    // vector<double> tanh;
    // for (int i = 0; i < tanh_encrypted.size(); i++) {
    //     vector<double> val;
    //     Plaintext val_plain;
    //     Ciphertext val_encrypted;
    //     val_encrypted = tanh_encrypted[i];
    //     decryptor.decrypt(val_encrypted, val_plain);
    //     encoder.decode(val_plain, val);
    //     tanh.push_back(val[0]);
    // }

    // print_vector(tanh, 500);
    
    ////////////////////////////////////////////////////////////////////////////
    // Decrypt and output data
    ////////////////////////////////////////////////////////////////////////////
    vector<double> wma12;
    for (int i = 0; i < wma12_encrypted_sliced.size(); i++) {
        vector<double> val;
        Plaintext val_plain;
        Ciphertext val_encrypted;
        val_encrypted = wma12_encrypted_sliced[i];
        decryptor.decrypt(val_encrypted, val_plain);
        encoder.decode(val_plain, val);
        wma12.push_back(val[0]);
    }

    vector<double> wma26;
    for (int i = 0; i < wma26_encrypted.size(); i++) {
        vector<double> val;
        Plaintext val_plain;
        Ciphertext val_encrypted;
        val_encrypted = wma26_encrypted[i];
        decryptor.decrypt(val_encrypted, val_plain);
        encoder.decode(val_plain, val);
        wma26.push_back(val[0]);
    }

    vector<double> wma_diff;
    for (int i = 0; i < wma_diff_sliced.size(); i++) {
        vector<double> val;
        Plaintext val_plain;
        Ciphertext val_encrypted;
        val_encrypted = wma_diff_sliced[i];
        decryptor.decrypt(val_encrypted, val_plain);
        encoder.decode(val_plain, val);
        wma_diff.push_back(val[0]);
    }

    vector<double> wma9;
    for (int i = 0; i < wma9_encrypted.size(); i++) {
        vector<double> val;
        Plaintext val_plain;
        Ciphertext val_encrypted;
        val_encrypted = wma9_encrypted[i];
        decryptor.decrypt(val_encrypted, val_plain);
        encoder.decode(val_plain, val);
        wma9.push_back(val[0]);
    }

    vector<double> macd;
    for (int i = 0; i < macd_encrypted.size(); i++) {
        vector<double> val;
        Plaintext val_plain;
        Ciphertext val_encrypted;
        val_encrypted = macd_encrypted[i];
        decryptor.decrypt(val_encrypted, val_plain);
        encoder.decode(val_plain, val);
        macd.push_back(val[0]);
    }
    
    // write data to csv for plotting
    ofstream output_file1("macd_seal.csv");
    ostream_iterator<double> output_iterator1(output_file1, "\n");
    copy(macd.begin(), macd.end(), output_iterator1);

    ofstream output_file2("wma12_seal.csv");
    ostream_iterator<double> output_iterator2(output_file2, "\n");
    copy(wma12.begin(), wma12.end(), output_iterator2);

    ofstream output_file3("wma26_seal.csv");
    ostream_iterator<double> output_iterator3(output_file3, "\n");
    copy(wma26.begin(), wma26.end(), output_iterator3);

    ofstream output_file4("wma_diff_seal.csv");
    ostream_iterator<double> output_iterator4(output_file4, "\n");
    copy(wma_diff.begin(), wma_diff.end(), output_iterator4);

    ofstream output_file5("wma9_seal.csv");
    ostream_iterator<double> output_iterator5(output_file5, "\n");
    copy(wma9.begin(), wma9.end(), output_iterator5);

    ofstream output_file6("decisions_macd.csv");
    ostream_iterator<double> output_iterator6(output_file6, "\n");
    copy(decisions.begin(), decisions.end(), output_iterator6);

    cout << "Data output done" << endl;

    return 0;
}
