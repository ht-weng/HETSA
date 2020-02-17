addpath('./Linear/Models')
addpath('./Linear/CorrelationFunctions')
addpath('./Linear/AutoregressiveModel')


%% Generate innovations
rng default
sampleSize = 1000;
inputSignal = (randn(1,sampleSize));
%% Generate timeseries using Moving Average Model
ma_coefs = [0.8 -0.1];
% Compare with Matlab standard function
rng default
Mdl = regARIMA('AR',{0.0},'MA', ma_coefs, 'Variance', 1, 'Intercept', 0);
ma1 = simulate(Mdl, sampleSize,'NumPaths', 1);
ma2 = maModel(inputSignal, ma_coefs);

% Plot both MA timeseries
figure, hold on,
plot(ma1, 'linewidth', 2);  
plot(ma2, 'linewidth', 2, 'Linestyle', ':');
legend({'MA1','MA2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Moving Average Model')

%% Generate timeseries using Autoregressive Model
ar_coefs = [1 -0.4 .2];
% Compare with Matlab standard function
rng default
Mdl = regARIMA('AR',ar_coefs,'MA', {0.0}, 'Variance', 1, 'Intercept', 0);
ar2 = simulate(Mdl, sampleSize,'NumPaths', 1);
ar3 = arModel(inputSignal, ar_coefs);

% Plot both AR timeseries
figure, hold on,
plot(ar2, 'linewidth', 2);  
plot(ar3, 'linewidth', 2, 'Linestyle', ':');
legend({'AR1','AR2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Autoregressive Model')

%% Generate timeseries using Autoregressive Moving Average Model
rng default
Mdl = regARIMA('AR',ar_coefs,'MA', ma_coefs, 'Variance', 1, 'Intercept', 0);
arma1 = simulate(Mdl, sampleSize,'NumPaths', 1);
arma2 = armaModel(inputSignal, ma_coefs, ar_coefs);

% Plot both AR timeseries
figure, hold on,
plot(arma1, 'linewidth', 2);  
plot(arma2, 'linewidth', 2, 'Linestyle', ':');
legend({'AR1','AR2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('ARMA Model')

% (2) Calculate covariance function
acf1 = autoCorrelation(arma1');
acf2 = autocorr(arma1, length(arma1) - 1);
disp(['Error in ACF estiomation:', num2str(norm(acf1 - acf2)) ]);
plotCorrelationFunction(acf1, 20);
plotCorrelationFunction(acf2, 20);


pacf1 = partialAutoCorrelation(ar3');
pacf2 = parcorr(ar3, length(ar3) - 1);
disp(['Error in PACF estiomation:', num2str(norm(pacf1(1:100) - pacf2(1:100)))]);
plotCorrelationFunction(pacf1, 100);
plotCorrelationFunction(pacf2, 100);