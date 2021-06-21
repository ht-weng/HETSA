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


// Recursive PID controller
// g only contains the G(t) and G(t-1)
inline tuple<Ciphertext, Ciphertext> pid_controller_recursive(vector<Ciphertext>& g, Ciphertext g_target, CKKSEncoder &encoder,
        Evaluator &evaluator, Decryptor &decryptor, RelinKeys &relin_keys, double scale, Ciphertext prev_sum) {

    // Encrypt parameters
    Plaintext u0, Kp, Ki_tor_i, Kd, neg_one;
    encoder.encode((16.67), scale, u0);
    encoder.encode((0.5), scale, Kp);
    encoder.encode((0.01), scale, Ki_tor_i);
    encoder.encode((0.05), scale, Kd);
    encoder.encode((-1), scale, neg_one);

    Ciphertext result_encrypted, gt_1, gt_2, first_term, second_term, third_term, cur_sum, dif;

    // G(t)
    gt_1 = g[1];
    // G(t-1)
    gt_2 = g[0];

    // e(t) = G(t) - g_target
    // e(t-1) = G(t-1) - g_target
    Ciphertext et_1, et_2;
    evaluator.sub(gt_1, g_target, et_1);
    evaluator.sub(gt_2, g_target, et_2);

    // Kp * e(t)
    evaluator.multiply_plain(et_1, Kp, first_term);
    evaluator.rescale_to_next_inplace(first_term);
    first_term.scale() = scale;

    // Ki * sum(e) / tor_i
    evaluator.mod_switch_to_inplace(prev_sum, et_1.parms_id());
    evaluator.add(et_1, prev_sum, cur_sum);
    evaluator.multiply_plain(cur_sum, Ki_tor_i, second_term);
    evaluator.rescale_to_next_inplace(second_term);
    second_term.scale() = scale;

    // Kd * (e(t) - e(t-1))
    // evaluator.multiply_plain(et_2, neg_one, third_term);
    // evaluator.rescale_to_next_inplace(third_term);
    // third_term.scale() = scale;
    // evaluator.mod_switch_to_inplace(et_1, third_term.parms_id());
    // evaluator.add_inplace(third_term, et_1);
    evaluator.sub(et_1, et_2, third_term);
    // evaluator.mod_switch_to_inplace(Kd, third_term.parms_id());
    evaluator.multiply_plain_inplace(third_term, Kd);
    evaluator.rescale_to_next_inplace(third_term);
    third_term.scale() = scale;

    // Add all
    evaluator.add(first_term, second_term, result_encrypted);
    // evaluator.mod_switch_to_inplace(result_encrypted, third_term.parms_id());
    evaluator.add_inplace(result_encrypted, third_term);
    evaluator.mod_switch_to_inplace(u0, result_encrypted.parms_id());
    evaluator.add_plain_inplace(result_encrypted, u0);

    // Return PID signal and current sum of errors
    return make_tuple(result_encrypted, cur_sum);
}


// Recursive PID controller
// g only contains the G(t) and G(t-1)
inline tuple<double, double> pid_controller_recursive_plain(vector<double>& g, double g_target, double prev_sum) {

    // Parameters
    double u0, Kp, Ki_tor_i, Kd, gt_1, gt_2, et_1, et_2, ut, cur_sum;
    u0 = 16.67;
    Kp = 0.5;
    Ki_tor_i = 0.01;
    Kd = 0.05;

    // G(t)
    gt_1 = g[1];
    // G(t-1)
    gt_2 = g[0];

    // e(t) = G(t) - g_target
    // e(t-1) = G(t-1) - g_target
    et_1 = gt_1 - g_target;
    et_2 = gt_2 - g_target;

    cur_sum = et_1 + prev_sum;

    ut = u0 + Kp * et_1 + Ki_tor_i * cur_sum + Kd * (et_1 - et_2);

    if (ut < 0) {
        ut = 0;
    }

    // Return ut and current sum of errors
    return make_tuple(ut, cur_sum);
}


// Bergman Minimal Model
inline vector<double> bergman(vector<vector<double>>& x, double U, double D) {
    // The last glucose signal G
    double G = x[x.size()-1][0];
    // The last insulin remote compartmen signal X
    double X = x[x.size()-1][1];
    // The last insulin signal I
    double I = x[x.size()-1][2];

    // Parameters
    double G0 = 4.5;
    double X0 = 15.0;
    double I0 = 15.0;
    // For T1D, Some papers P1 = 0 to T1D
    double P1 = 0.028735;
    double P2 = 0.028344;
    double P3 = 5.035e-05;
    double VI = 12.0;
    double n = 0.09259259;
    // Minimal Model
    double Gdt = -P1 * (G - G0) - (X - X0) * G + D;
    double Xdt = -P2 * (X - X0) + P3 * (I - I0);
    double Idt = -n * I + U/VI;

    vector<double> dx_dt{Gdt, Xdt, Idt};
    return dx_dt;
}

// Create a meal profile to simulate one afternoon of diabetes patient
inline double meal_profile(int t) {
    double m;
    if (t < 100) {
        m = 0.0;
    } else if (t >= 100 & t < 9000){
        // Lunch
        m = 2.0 * exp(double(-0.001 * double(t-100)));
    } else if (t >= 9000 & t < 16000) {
        // Snack
        m = 1.0 * exp(double(-0.001 * double(t-9000)));
    } else if (t >= 16000 & t < 21600) {
        // Dinner
        m = 2.0 * exp(double(-0.001 * double(t-16000)));
    } else {
        m = 0.0;
    }
    return m;
}

int main() {
    // SEAL settings
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    double scale = pow(2.0, 40);
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

    // Apply recursive PID on Bergman Minimal Model
    // Only the last two errors and the previous sum are passed into the algorithm
    cout << "Applying recursive PID controller" << endl;
    auto recpid_start = high_resolution_clock::now();

    // Total time (seconds)
    int time_limit = 21600;

    // Set up initial conditions
    vector<double> x0{4.5, 15.0, 15.0};
    vector<double> x1{4.5, 15.0, 15.0};
    vector<vector<double>> x{x0, x1};

    // Set point of blood glucose level
    double g_target = 6.0;
    Plaintext g_target_plain;
    Ciphertext g_target_encrypted;
    encoder.encode(g_target, scale, g_target_plain);
    encryptor.encrypt(g_target_plain, g_target_encrypted);

    // Meal profiles
    vector<double> m;
    for (int i = 0; i < time_limit; i++) {
        m.push_back(meal_profile(i));
    }

    // Initialise previous sum
    double zero = 0.0;
    Plaintext zero_plain;
    Ciphertext zero_encrypted;
    encoder.encode(zero, scale, zero_plain);
    encryptor.encrypt(zero_plain, zero_encrypted);
    Ciphertext prev_sum_encrypted;
    prev_sum_encrypted = zero_encrypted;


    // Simulation
    for (int t = 0; t < time_limit; t++) {
        vector<double> dx_dt;
        tuple<Ciphertext, Ciphertext> tup;
        Ciphertext ut_encrypted;
        Plaintext ut_plain;
        vector<double> ut;
        
        // Get G(t) and G(t-1)
        double gt1 = x[x.size()-1][0];
        double gt2 = x[x.size()-2][0];

        cout << "Time point (second): " << t << endl;
        cout << "Current glucose level: " << gt1 << endl;

        // Encrypt G(t) and G(t-1) and store them into a vector
        Plaintext gt1_plain, gt2_plain;
        encoder.encode(gt1, scale, gt1_plain);
        encoder.encode(gt2, scale, gt2_plain);
        Ciphertext gt1_encrypted, gt2_encrypted;
        encryptor.encrypt(gt1_plain, gt1_encrypted);
        encryptor.encrypt(gt2_plain, gt2_encrypted);
        vector<Ciphertext> g_encrypted{gt2_encrypted, gt1_encrypted};

        // Apply PID on encrypted blood glucose level
        tup = pid_controller_recursive(g_encrypted, g_target_encrypted, encoder, evaluator, decryptor, 
            relin_keys, scale, prev_sum_encrypted);
        
        // Decrypt ut
        ut_encrypted = get<0>(tup);
        decryptor.decrypt(ut_encrypted, ut_plain);
        encoder.decode(ut_plain, ut);

        // Update previous sum
        prev_sum_encrypted = get<1>(tup);

        // Bergman Minimal Model simulation
        dx_dt = bergman(x, ut[0], m[t]);
        
        // x(t) = x(t-1) + dx/dt
        vector<double> xt;
        for (int i = 0; i < 3; i++) {
            xt.push_back(x[x.size()-1][i]+dx_dt[i]);
        }
        
        // Append xt to x
        x.push_back(xt);
    }

    // double prev_sum = 0.0;
    // vector<double> dx_dt;
    // tuple<double, double> tup;
    // double ut;

    // // Simulation
    // for (int t = 0; t < time_limit; t++) {

    //     // Get G(t) and G(t-1)
    //     double gt1 = x[x.size()-1][0];
    //     double gt2 = x[x.size()-2][0];

    //     // cout << "gt: " << gt << endl; 
    //     // cout << "gt1: " << gt1 << endl;

    //     vector<double> g{gt2, gt1};

    //     // Apply PID on encrypted blood glucose level
    //     tup = pid_controller_recursive_plain(g, g_target, prev_sum);
        
    //     ut = get<0>(tup);

    //     // Update previous sum
    //     prev_sum = get<1>(tup);

    //     // Bergman Minimal Model simulation
    //     dx_dt = bergman(x, ut, m[t]);
        
    //     // x(t) = x(t-1) + dx/dt
    //     vector<double> xt;
    //     for (int i = 0; i < 3; i++) {
    //         xt.push_back((x[x.size()-1][i]+dx_dt[i]));
    //     }
        
    //     // Append xt to x
    //     x.push_back(xt);
    // }

    auto recpid_stop = high_resolution_clock::now();
    auto recpid_duration = duration_cast<seconds>(recpid_stop - recpid_start);
    cout << "Recursive PID Duration:  " << recpid_duration.count() << " seconds" << endl;
    cout << endl;

    cout << "Output data" << endl;
    cout << endl;

    vector<double> G, X, I;
    G = get_column(x, 0);
    X = get_column(x, 1);
    I = get_column(x, 2);
    exportData(G, "../../data/pid_bergman_recursive_local_G.csv");
    exportData(X, "../../data/pid_bergman_recursive_local_X.csv");
    exportData(I, "../../data/pid_bergman_recursive_local_I.csv");
    exportData(m, "../../data/meal_profile.csv");
    cout << endl;
}
