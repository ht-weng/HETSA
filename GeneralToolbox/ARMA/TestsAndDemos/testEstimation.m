addpath('./Models')
addpath('./CorrelationFunctions')
addpath('./AutoregressiveModel')

%% Generate innovations
rng default
sampleSize = 1000;
inputSignal = (randn(1,sampleSize));

% Test AR estimation using Durbin Levinson
ar_coefs = [1 - 0.4 .2];
rng default
ar_process = arModel(inputSignal, ar_coefs);

% Estimate parameters
[ar_coefs_est, error, pacf, significance ]= estimateARbyDurbinLevinson(ar_process, length(ar_coefs));

ar1 = arModel(inputSignal, ar_coefs);
ar2 = arModel(inputSignal, ar_coefs_est);
% Plot both AR timeseries
figure, hold on,
plot(ar1, 'linewidth', 2);  
plot(ar2, 'linewidth', 2, 'Linestyle', ':');
legend({'AR1','AR2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Autoregressive Model')

% Test AR estimation using Burg
[ar_coefs_est, error, pacf, significance ]= estimateARbyBurgs(ar_process, length(ar_coefs));

ar1 = arModel(inputSignal, ar_coefs);
ar2 = arModel(inputSignal, ar_coefs_est);
% Plot both AR timeseries
figure, hold on,
plot(ar1, 'linewidth', 2);  
plot(ar2, 'linewidth', 2, 'Linestyle', ':');
legend({'AR1','AR2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Autoregressive Model')

% Test MA estimation using Burg
ma_coefs = [0.8 -0.1];
rng default
ma_process = maModel(inputSignal, ma_coefs);
[ma_coefs_est, error] = estimateMAbyInnov(ma_process, length(ma_coefs));

ma1 = maModel(inputSignal, ma_coefs);
ma2 = maModel(inputSignal, ma_coefs_est);

% Plot both MA timeseries
figure, hold on,
plot(ma1, 'linewidth', 2);  
plot(ma2, 'linewidth', 2, 'Linestyle', ':');
legend({'MA1','MA2'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Moving Average Model')