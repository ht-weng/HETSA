# Homomorphically Encrypted Time Series Analysis (HETSA)

## Statement 
The primary focus of the project is to develop and implement tools for storing, representing, and analysing encrypted time-series data. Furthermore, the proposed representation and bootstrapping algorithm will allow for real-time analysis and decision making. The expected impact of this research project is twofold. Firstly, the proposed framework will lay a foundation for further research and development of methods of time series analysis. The developed tools will serve as building blocks for more sophisticated and advanced machine learning algorithms such as recurrent neural networks. Secondly, the proposed framework implemented as (HETSA) open-source library will be ready to use and applied for domain-specific time-series data. This in turn will ignite further research as well open up new application areas.

## Current problems and limmitations

The idea of homomorphic encryption was proposed in the 70s as a solution to privacy preserving computation. Modern cryptosystems are capable of performing arbitrary computation on ciphertexts, facilitating the implementation of various machine learning algorithms in the homomorphic encryption framework. This is an active area of research, and therefore a multitude of encryption schemes have been proposed and implemented as open-source libraries.   Nevertheless, there are a number of limitations associated with existing libraries:

1. There is a multitude of open source libraries that often implement equivalent homomorphic encryption schemes. This variety of libraries complicates the decision of choosing the right one.
2. These libraries are relatively low-level and implement elementary operations of multiplication and addition. Hence, making the implementation of data analysis algorithms arduous.
3. Scalars, vectors and matrices are mapped to sets of coefficients that make the union operation over encrypted datasets either impossible or tedious. This in turn obfuscate the process of real-time time series analysis.

## Developers

Haotian Weng, Hanna Suominen, Artem Lenskiy


# Applications
To demonstrate the capabilities of the proposed framework we present a number of use cases.

## 1. Moving Average Convergence Divergence with SEAL

In this application, we demonstrate an analogous to Numer.ai system implemented with the HTSA library. In our scenario a set of traders receive encrypted price quotes, analyse the prices in real time, and draw a trading decision (if any). The trading decision in encrypted form is sent back to the server where it is decrypted and facilitated.

## 2. Diabetes diagnosis on encrypted heart rate variability

In this application, HETSA library is employed to detect diabetes through the analysis of encrypted heart rate variability (HRV). HRV is usually measured by a wearable device (e.g. smartwatch) and is transmitted over the Internet in an encrypted form to a remote server for further processing. In order to make a diagnosis, the HETSA implements spectral analysis and linear classifiers. 

# Building

### Dependencies
Microsoft SEAL library
cmake
gcc/g++
### Getting started
````
cmake .
make
./main
````
