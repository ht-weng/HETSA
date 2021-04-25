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
#include <tuple>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <iterator>
#define PORT 8080
# include "seal/seal.h"

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

// Original PID controller
// error_hist consists of all previous errors
inline Ciphertext pid_controller(vector<Ciphertext>& error_hist, CKKSEncoder &encoder,
        Evaluator &evaluator, Decryptor &decryptor, RelinKeys &relin_keys, double scale) {
        
    // Encrypt parameters
    Plaintext u_steady, Kp, Ki_tor_i, Kd, neg_one;
    encoder.encode((16.67), scale, u_steady);
    encoder.encode((0.5), scale, Kp);
    encoder.encode((0.01), scale, Ki_tor_i);
    encoder.encode((0.05), scale, Kd);
    encoder.encode((-1), scale, neg_one);

    Ciphertext result_encrypted, last_signal, sec_last_signal, first_term, second_term, third_term, cur_sum, dif;

    // error_hist[-1]
    last_signal = error_hist[error_hist.size()-1];
    // error_hist[-2]
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

    // Return PID signal
    return result_encrypted;
}

// Recursive PID controller
// error_hist only consists of the last two errors
inline tuple<Ciphertext, Ciphertext> pid_controller_recursive(vector<Ciphertext>& error_hist, CKKSEncoder &encoder,
        Evaluator &evaluator, Decryptor &decryptor, RelinKeys &relin_keys, double scale, Ciphertext prev_sum) {
        
    // Encrypt parameters
    Plaintext u_steady, Kp, Ki_tor_i, Kd, neg_one;
    encoder.encode((16.67), scale, u_steady);
    encoder.encode((0.5), scale, Kp);
    encoder.encode((0.01), scale, Ki_tor_i);
    encoder.encode((0.05), scale, Kd);
    encoder.encode((-1), scale, neg_one);

    Ciphertext result_encrypted, last_signal, sec_last_signal, first_term, second_term, third_term, cur_sum, dif;

    // error_hist[-1]
    last_signal = error_hist[error_hist.size()-1];
    // error_hist[-2]
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

    // Return PID signal and current sum of errors
    return make_tuple(result_encrypted, cur_sum);
}


// Bergman Minmod
// Human Glucose - Insulin System Model by Bergman 1981
// A simple 3 equation model.

// States = [Plasma Glucose(G), 
//           Plasma Insulin Remote Compartment(X), 
//           Plasma Insulin(I)]
inline double bergman_minmod(double U) {
    // dummy bergman for testing
    return U;
}


int main() {
    // int SOCKET_SIZE = 2097152;
    // Socket settings
    int obj_socket = 0, reader;
    struct sockaddr_in serv_addr;
    // Set buffer size
    char buffer[2097152] = {0};
    if ((obj_socket = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        cout << "Socket creation error !" << endl;
        return -1;
    }
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(PORT);
    // Converting IPv4 and IPv6 addresses from text to binary form
    if (inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr) <= 0)
    {
        cout << "Invalid address ! This IP Address is not supported !" << endl;
        return -1;
    }
    if (connect(obj_socket, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        cout << "Connection Failed : Can't establish a connection over this socket!" << endl;
        return -1;
    }

    // Create stringstreams
    stringstream parms_stream;
    stringstream data_stream;
    stringstream sk_stream;

    // SEAL settings
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));
    double scale = pow(2.0, 40);
    SEALContext context(parms);
    print_parameters(context);
    cout << endl;

    // Send Parms
    auto parms_size = parms.save(parms_stream);
    cout << "Parms size: " << parms_size << " bytes" << endl;
    send(obj_socket, parms_stream.str().c_str(), parms_stream.str().length(), 0);
    cout << "Parms sent to Bergman" << endl;

    // Read relin keys from bergman
    reader = read(obj_socket, buffer, 2097152);
    data_stream << buffer;
    RelinKeys relin_keys;
    relin_keys.load(data_stream);
    cout << "Load relin key from bergman" << endl;

    // Reset data stream
    data_stream.str(string());

    // Set up encoder, encryptor, evaluator and decryptor
    CKKSEncoder encoder(context);
    Evaluator evaluator(context);

    // Recursive PID with Socket
    // Only the last two errors and the previous sum are passed into the algorithm
    // cout << "Applying recursive PID controller on Socket" << endl;
    // auto recpid_start = high_resolution_clock::now();

    vector<Ciphertext> result_rec_socket_encrypted;
    vector<Ciphertext> prev_sum_vec_socket;
    // Initialisation
    Ciphertext init_sum = error_hist[0];
    // Vector of the first two errors
    vector<Ciphertext> sliced_hist_rec = slice(error_hist, 0, 2);
    // Generate the first signal
    tuple<Ciphertext, Ciphertext> tup_socket = pid_controller_recursive(sliced_hist_rec, encoder, evaluator, decryptor, relin_keys, scale, init_sum);
    // Store the first signal and sum of the first two errors
    result_rec_socket_encrypted.push_back(get<0>(tup_socket));
    prev_sum_vec_socket.push_back(get<1>(tup_socket));



    // Recursively generate PID signals
    // The last two errors and the previous sum are passed into the algorithm
    for (int i = 3; i < error_hist.size(); i++) {

        // sliced_hist_rec = slice(error_hist, i-2, i);

        // // Read hist_rec from server
        // reader = read(obj_socket, &sliced_hist_rec, sizeof(sliced_hist_rec));

        Signal sliced_hist_rec;

        // Read the existing address book.
        fstream input(argv[1], ios::in | ios::binary);
        if (!input)
        {
            cout << argv[1] << ": File not found.  Creating a new file." << endl;
        }
        else if (!sliced_hist_rec.ParseFromIstream(&input))
        {
            cerr << "Failed to parse address book." << endl;
            return -1;
        }

        tup_socket = pid_controller_recursive(sliced_hist_rec, encoder, evaluator, decryptor, 
            relin_keys, scale, prev_sum_vec_socket[prev_sum_vec_socket.size()-1]);
        
        // // Decrypt U and feed it into Bergman Minmod
        // Ciphertext U_encrypted = get<0>(tup_socket);
        // Plaintext U_plaintext;
        // double U;
        // decryptor.decrypt(U_encrypted, U_plaintext);
        // encoder.decode(U_plaintext, U);
        // double dx_dt = bergman_minmod(U);

        // // Encrypt
        // Ciphertext X_encrypted;

        // // Send tup_socket to server
        // send(obj_socket, &tup_socket, sizeof(tup_socket), 0);

        // Write the new address book back to disk.
        fstream output(argv[1], ios::out | ios::trunc | ios::binary);
        if (!tup_socket.SerializeToOstream(&output))
        {
            cerr << "Failed to write address book." << endl;
            return -1;
        }

        // result_rec_socket_encrypted.push_back(get<0>(tup_socket));
        prev_sum_vec_socket.push_back(get<1>(tup_socket));
    }
    auto recpid_stop = high_resolution_clock::now();
    auto recpid_duration = duration_cast<seconds>(recpid_stop - recpid_start);
    cout << "Recursive PID on socket Duration:  " << recpid_duration.count() << " seconds" << endl;
    cout << endl;

    // cout << "Output data" << endl;
    // cout << endl;

    // vector<double> result_rec_socket_plaintext = decryptVec(result_rec_socket_encrypted, decryptor, encoder);
    // double rpid_error = pe(signal_input, result_rec_socket_plaintext);
    // cout << "Error between recursive PID and plaintext PID: " << rpid_error*100 << "%" << endl;
    // exportData(result_rec_socket_plaintext, "../data/seal_result_recursive_socket.csv");
    // cout << endl;
}
