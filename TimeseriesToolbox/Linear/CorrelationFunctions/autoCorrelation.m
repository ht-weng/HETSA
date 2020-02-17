function [gamma, significance] = autoCorrelation(x)
    if(size(x, 1) > size(x, 2))
         disp('Make sure that timeseries are stored as rows');
         gamma = [];
         significance = [];
         return;
    end
    % remove the mean value
    x = x - mean(x);
    lengthX = length(x);
    % augment X with zeros
    extX = [x, zeros(1,lengthX)];
    gamma = zeros(1, lengthX);
    %caclulate the covariance sum_t ( X(t)*X(t+h) ) / n
    for h = 1:lengthX
        gamma(h) = sum(x .* extX(h:h+lengthX - 1));
    end
    gamma = gamma / (lengthX - 1);
    gamma = gamma/gamma(1); % divide by variace
    significance  = (abs(gamma) > (1.96 / sqrt(lengthX)));
end