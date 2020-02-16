addpath('./Linear/Models')
addpath('./Linear/CorrelationFunctions')
addpath('./Linear/AutoregressiveModel')

%% Generate timeseries to be tested
t = linspace(0,2*pi,100);
rng default  %initialize random number generator
x = [1*t; 0*t] +...                                             % Linear trend
    [1; -1].*sin([1; 0.5].*t)...                                % Nonlinear trend
    + sin(5*t)...                                              % Periodic component
    + [0.25;0.1].*randn(2, size(t, 2))...                         % Gaussian noise
    + (rand(2, size(t, 2)) > 0.9).*[1;2].*randn(2, size(t, 2))   % Impuse noise
    + 5;                                                        % Bias

%% Apply Exponential Moving Average
% Prepare the coefficients
alpha = 0.35;
b = [alpha];
a = [1  -(1-alpha)];
% Visualize filter frequecny respone
[h,w] = freqz(b,a,128);
figure, plot(w/pi,abs(h))
xlabel 'Radian frequency (\omega/\pi)', ylabel Magnitude


y1 = filter(b, a, x, [], 2);
y2 = armaFilter(b, a, x);

figure, hold on;
plot(t,x)
plot(t, y1, '-.')
plot(t, y2, 'linewidth', 2)
legend({'x_1','x_2', 'y_1(x_1)', 'y_1(x_2)', 'y_2(x_1)', 'y_2(x_2)'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Exponential Moving Average')
%% Apply Recursive digital filter design
% Prepare the coefficients
[b, a] = yulewalk(5, [0 0.5 0.5 1], [1 1 0 0]);
% Visualize filter frequecny respone
[h,w] = freqz(b, a, 128);
% figure, plot(w/pi,abs(h),f,m,'--')
% xlabel 'Radian frequency (\omega/\pi)', ylabel Magnitude
% legend('Yule-Walker','Ideal'), legend boxoff

y1 = filter(b,a,x);
y2 = armaFilter(b, a, x);

figure, hold on
plot(t, x)
plot(t, y1, '-.')
plot(t, y2, 'linewidth', 2)
%plot([7.25929,7.40171,7.19729,7.66486,7.59457,7.22643,7.42129,7.76786,8.374,8.52129,7.39057,8.80671,8.92043,8.11357,8.85771,7.46714,8.16714,9.49757,9.29314,9.88429])
legend({'x_1','x_2', 'y_1(x_1)', 'y_1(x_2)', 'y_2(x_1)', 'y_2(x_2)'}, 'Orientation', 'vertical', 'Location', 'NorthWest');
title('Recursive digital filter');
