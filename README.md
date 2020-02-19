# Homomorphically Encrypted Time Series Analysis  (HETSA)  

## Statement

The aim of this project is to develop a framework for analysing homomorphically encrypted time series (HETSA). Namely, we will implement a set of building blocks for time series analysis and showcase them on the following two applications: First, to demonstrate the capabilities of real-time processing, we will implement a set of trading strategies that make trading decisions on encrypted market data. Second, to demonstrate the developed signal analysis tools, we will analyse encrypted heart rate variability of patients with diabetes that consequently will support clinical judgement and decision-making. Third, to show the versatility, we will implement a set of general time-series analysis tools in homomorphic encryption domain. Furthermore, we will develop a bootstrapping technique that provides the means for real time data processing. We believe that focusing on these areas, where issues surrounding privacy play a crucial role, will not only showcase the potential of the proposed framework, but also ignite further research and development.  

## Current problems and limmitations

The idea of homomorphic encryption was proposed in the 1970s as a solution to privacy-preserving computation. Modern cryptosystems are capable of performing arbitrary computation on ciphertexts, facilitating the implementation of various machine learning algorithms in the homomorphic encryption context. This is an active area of research, and therefore a multitude of encryption schemes have been proposed and implemented as open-source libraries.  

Nevertheless, existing libraries have a number of limitations: First, an array of open source libraries that implement equivalent homomorphic encryption schemes are available to choose from, and thereby, complicate the decision to choose the most relevant one. Second, these libraries are relatively low-level and implement elementary operations such as multiplication and addition, hence, making the implementation of data analysis algorithms arduous. Third, scalars, vectors, and matrices are mapped to sets of coefficients that make the union and intersection operations over encrypted datasets either impossible or tedious. Thus, making the process of real-time time series analysis cumbersome. Last but not least, as the state-of-art bootstrapping method in CKKS scheme is inefficient and impractical, at the moment the multiplication level of our algorithms is restricted by the parameters of the CKKS scheme and thus it is difficult to implement recursive algorithms which is a foundation of timeseries analysis.  

## Developers

Haotian Weng, Hanna Suominen, Artem Lenskiy

## Applications

To demonstrate the capabilities of the proposed framework we present a number of use cases.

GeneralToolbox
    ├── ARMA - To be implemented in homomorphic encryption context
    │
    ├── FinancialApplications
    │   ├── macd-HEAAN
    │   │
    │   ├── macd-Plaintext
    │   │
    │   └── macd-SEAL
    │
    └── MedicalApplications - TODO

### General Toolbox - Autoregressive Moving Average Model

ARMA model, which is a popular time-series analysis tool, is the combination of Autoregressive (AR) model and Moving Average (MA) model. AR model is a linear regression model which predicts future values based on the history data and works well when there is correlation between the values in the time-series. MA model generates a value linearly dependent on the current and past values, which indicates the trend or momentum of the time-series data. Combining the advantages of AR modle with the ones of MA model, ARMA model usually requires fewer parameters than AR or MA model alone.  With the implementation of ARMA model in homomorphic encryption context, time-series analysis on various encrypted data becomes possible, which is reason why it is a general-purpose tool for HETSA.  

### Financial Applications - Moving Average Convergence Divergence

Algorithmic trading has proliferated the area of quantitative finance for already a number of decades. The decisions are made in an automatic manner using the data provided by brokerage firms and exchanges. There is an emerging intermediate layer of financial players that is placed in between a broker and algorithmic traders. The role of these players is to aggregate market decisions from the algorithmic traders and send a final market order to a broker. In return the quants receive incentives proportional to the correctness of their predictions. In such a setup, the intermediate player - an aggregator does not provide the market data in plaintext but encrypts it. Encrypting market data prevents quants from trading on their own, as well as keeps expensive financial data private. In this use case scenarion we implement a MACD-based trend-following strategy using the methods of homomorphic encryption.

### Medical Applications - Diabetes diagnosis on encrypted heart rate variability - TODO

In this application, the HETSA library is employed to detect diabetes through the analysis of the encrypted heart rate variability (HRV). HRV is usually measured by a wearable device (e.g. smartwatch) and is transmitted over the Internet in an encrypted form to a remote server for further processing. In order to make a diagnosis, the HETSA implements spectral analysis and linear classifiers.  

## Getting started  

### General Toolbox

#### ARMA

Run the test cases in ARMA/TestsAndDemos in MATLAB.

### Financial Application

Developed and tested on macOS Catalina.  
Please note that the minimum RAM required is 32GB as Ciphertexts are memory-consuming.  

#### Dependencies

[Microsoft SEAL Library](https://github.com/Microsoft/SEAL)  
[HEAAN Library](https://github.com/snucrypto/HEAAN)  
cmake  
gcc/g++  

#### MACD with SEAL  

First install [Microsoft SEAL Library](https://github.com/Microsoft/SEAL)  

```bash
cd macd-SEAL/
cmake .
make
./macd
```

#### MACD with HEAAN

HEAAN Library is included in macd-HEAAN.  

```bash
cd macd-HEAAN/run
make
./MACD
```

#### MACD in Plaintext

Run macd.ipynb in Jupyter Notebook to generate MACD signals and trading decisions based on plaintext data.  

#### Visualisation

Run output-comparison.ipynb in Jupyter Notebook to visualise results.  
