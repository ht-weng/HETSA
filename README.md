# Analysis of Homomorphically Encrypted Time Series  (HETSA)

## Statement 
The aim of this project is to develop a framework for analysing homomorphically encrypted time series (HETSA). Namely, we will develop a method of packing encrypted data points into a vector and a bootstrapping technique that provide the means for real time data processing. We will further implement a set of building blocks for time series analysis and showcase them on the following two applications: First, to demonstrate the capabilities of real-time processing, we will implement a set of trading strategies that make trading decisions on encrypted market data. Second, to demonstrate the developed signal analysis tools, we will analyse encrypted heart rate variability of patients with diabetes that consequently will support clinical judgement and decision-making. We believe that focusing on these two major areas, where issues surrounding privacy play a crucial role, will not only showcase the potential of the proposed framework, but also ignite further research and development. 

## Current problems and limmitations

The idea of homomorphic encryption was proposed in the 1970s as a solution to privacy-preserving computation. Modern cryptosystems are capable of performing arbitrary computation on ciphertexts, facilitating the implementation of various machine learning algorithms in the homomorphic encryption context. This is an active area of research, and therefore a multitude of encryption schemes have been proposed and implemented as open-source libraries.  
Nevertheless, existing libraries have a number of limitations: First, an array of open source libraries that implement equivalent homomorphic encryption schemes are available to choose from, and thereby, complicate the decision to choose the most relevant one. Second, these libraries are relatively low-level and implement elementary operations of multiplication and addition, hence, making the implementation of data analysis algorithms arduous. Third, scalars, vectors, and matrices are mapped to sets of coefficients that make the union and intersection operations over encrypted datasets either impossible or tedious. Thus, making the process of real-time time series analysis cumbersome.



## Developers

Haotian Weng, Hanna Suominen, Artem Lenskiy


# Applications
To demonstrate the capabilities of the proposed framework we present a number of use cases.

## 1. Moving Average Convergence Divergence with SEAL

In this application, we demonstrate an analogous to Numer.ai system implemented with the HETSA library. In our scenario a set of traders receive encrypted price quotes, analyse the prices in real time, and draw a trading decision (if any). The trading decision in encrypted form is sent back to the server where it is decrypted and facilitated.

## 2. Diabetes diagnosis on encrypted heart rate variability

In this application, the HETSA library is employed to detect diabetes through the analysis of the encrypted heart rate variability (HRV). HRV is usually measured by a wearable device (e.g. smartwatch) and is transmitted over the Internet in an encrypted form to a remote server for further processing. In order to make a diagnosis, the HETSA implements spectral analysis and linear classifiers. 

# Building

### Dependencies
[Microsoft SEAL Library](https://github.com/Microsoft/SEAL)   
[HEAAN Library](https://github.com/snucrypto/HEAAN)  
cmake  
gcc/g++  
### Getting started  
Developed and tested on macOS Catalina.  
Please note that the minimum RAM required is 16GB as Ciphertexts are memory-consuming.  
#### MACD with SEAL:  
First install [Microsoft SEAL Library](https://github.com/Microsoft/SEAL)  
````
cd macd-SEAL/
cmake .
make
./macd
````
Run macd-seal.ipynb to visualise results.  
#### MACD with HEAAN:  
HEAAN Library is included in macd-HEAAN.  
````
cd macd-HEAAN/run
make
./MACD
````
Run macd-heaan.ipynb to visualise results.  
